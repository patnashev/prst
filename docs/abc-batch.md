# ABC batch mode — deep dive

`-batch` runs PRST over a *list* of candidates instead of one, driving the same `Run::create → run->run()` machinery once per entry, with a single shared `Options` config, cross-candidate progress/resume, and stop conditions. The list comes from a file (or `stdin`) in one of five shapes: **raw** (one expression per line), **ABC**, **ABCD**, **ABC2**, or **NewPGen** — the output formats of the common sieves (srsieve2, gcwsieve, mfsieve, NewPGen). Two cleanly separated pieces:

- `abc_parser.{h,cpp}` — pure parsing: turn a batch file into a `CandidateSource`, an indexable sequence of candidate expression strings. No primality logic.
- `batch.{h,cpp}` — the driver: `batch_main` parses options, opens the source, and loops, calling the same `InputNum`/`Run` path the single-candidate `main()` uses.

`batch_main` is an **alternate entry point**, not a mode inside `main()`: `prst.cpp:214` does `.check_code("-batch", [&]{ exit(batch_main(argc, argv)); })` — the same `exit()`-into-a-sibling pattern as `-boinc`/`-net`. It re-uses the test classes; it doesn't replace them.

Source files:
- `src/abc_parser.h` (`FileFormat`, `Candidate`, `CandidateSource`, `parse_batch_file`, the limit constants)
- `src/abc_parser.cpp` (format detection, `ABCTemplate`, the two `CandidateSource` impls, the per-format parsers)
- `src/batch.cpp` (`batch_main` — option parsing + the driver loop)
- `src/prst.cpp:214` (the `-batch` dispatch)

Prereqs / companions: `inputnum-parsing.md` (each candidate string is fed to `InputNum::parse`; the small-number path calls `factorize_small`), `run-hierarchy.md` (`Run::create` + `run->run()` per candidate), `logging-and-progress.md` (the batch keeps a *second* `Logging` for cross-candidate progress and the `cur`/`primes`/`composites` resume params). The format module was split out in commit `2722d20`; the `-stop` conditions were added in `020855b`.

## 1. Data model

The parser's job is to produce a `CandidateSource` — a flat, indexable list of `Candidate`s:

```cpp
struct Candidate
{
    std::string expression;   // e.g. "5*2^100+1", fed verbatim to InputNum::parse
    std::string k_value;      // the "k" for per-k tracking; empty for non-ABC formats
};

class CandidateSource
{
public:
    virtual size_t size() const = 0;
    virtual bool get(size_t index, Candidate& out) const = 0;   // false on bad index
    virtual bool is_abc() const = 0;                             // has k-values?
};
```

Two implementations (`abc_parser.cpp`):

- **`VectorCandidateSource`** (`:435-464`) — materialized: holds `vector<string>` of expressions and k-values. Used for **ABC, ABCD, NewPGen, and raw**. `get` is a bounds-checked array read.
- **`ABC2CandidateSource`** (`:472-535`) — **lazy in the product, not the factors**: it synthesizes each candidate on demand from a flat index via mixed-radix decomposition (`_strides`, innermost = last variable). The per-variable value *lists* are materialized (a `primes from 1 to 10^8` range becomes a `vector` of ~5.7M `int64`s), but the **Cartesian product** is not — memory is O(sum of range sizes), not O(product), so a product of millions of candidates costs only the per-variable lists.

Supporting structs:

- **`ABCTemplate`** (`:154-286`) — a parsed `$a*2^$b+1`-style template: alternating `literals` and `var_indices` (`$a`…`$d` → `0`…`3`), plus the **k-variable** (the `$var` or fixed text before the first `*`; `"1"` if there's no `*`). `expand(data_line, …)` substitutes whitespace-separated values into the template **verbatim** — including any sign character on a value token (see §5 on the `+$c`-vs-`$c` convention).
- **`ABC2VarRange`** (`:292-296`) — one variable's `{var_index, values}`. `parse_abc2_var_line` (`:298`) fills `values` from an explicit list, an arithmetic range, or a prime sieve.
- **`FileFormat`** enum (`abc_parser.h:39`) + `detect_format(first_line)` (`abc_parser.cpp:86-106`) — case-insensitive prefix match on `"ABCD "` / `"ABC2 "` / `"ABC "` (note the **trailing space** — `"ABCD"` with no space is `UNKNOWN`), then a NewPGen header sniff: a line that parses as ≥ 4 colon-separated fields (`<sievelimit>:<char>:<chainlen>:<base>[:<mask>]`) is `FORMAT_NEWPGEN`. Anything else is `UNKNOWN` → raw.

Hard limits (`abc_parser.h:27-33`), all guard against a malformed/hostile sieve file blowing up memory:

| Constant | Value | Guards |
|---|---|---|
| `MAX_SIEVE_RANGE` | 100,000,000 | the span of an ABC2 `primes from X to Y` sieve |
| `MAX_ABC2_CANDIDATES` | 100,000,000 | the total ABC2 Cartesian product |
| `MAX_ABC2_IN_LIST_VALUES` | 1,000 | values in an ABC2 `in { … }` list (excess is **truncated**, with a warning) |

## 2. `batch_main` — the central function

`batch.cpp:31-388`. After option parsing (a `Config` DSL subset of `main()`'s — `:48-138`) and the usage dump when no batch name is given, the driver is one big loop. The skeleton, annotated:

```cpp
source = parse_batch_file(batch_name, logging_batch);        // unless batch_name == "stdin"
size_t total = source ? source->size() : 0;

// resume: a SECOND progress file, keyed by batch name + base/order suffix
File batch_progress(batch_name + filename_suffix + ".param", 0);
logging_batch.file_progress(&batch_progress);
int cur = logging_batch.progress().param_int("cur");          // resume index
int primes     = logging_batch.progress().param_int("primes");
int composites = logging_batch.progress().param_int("composites");
std::map<std::string, bool> k_prime_found;

bool success = false;
for (; batch_name == "stdin" || cur < (int)total; cur++)
{
    logging_batch.report_param("cur", cur);
    logging_batch.progress().update(total > 0 ? cur/(double)total : 0, 0);
    logging_batch.progress_save();

    if (success && stop_prime)        { Task::abort(); break; }      // stop AFTER a prime
    if (stop_composites > 0 && composites >= stop_composites) { … break; }

    // 1. fetch the next expression (stdin getline, or source->get(cur))
    // 2. stop_k_prime: skip if a prime was already found for this k_value
    InputNum input;
    if (!input.parse(expression)) { /* print error; stop_error? abort : continue */ }

    if (input.bitlen() <= 40) { /* factorize_small → prime/not-prime; continue */ }
    else if (trial_division) { /* factorize_small; if factor → not-prime; continue */ }

    // per-candidate files: prst_<fingerprint><suffix>.{param,ckpt,rcpt}
    File file_progress(filename_prefix + filename_suffix + ".param", fingerprint);
    logging.file_progress(&file_progress);                            // load THIS candidate's params

    std::unique_ptr<Run> run(Run::create(input, options, logging));   // per-candidate logging
    if (!run) continue;

    GWState gwstate;                                                  // fresh GWState per candidate
    options.configure(gwstate);                                       // thread/cpu/fft knobs from the shared Options
    if (gwstate.next_fft_count < logging.progress().param_int("next_fft"))
        gwstate.next_fft_count = logging.progress().param_int("next_fft");  // restore a persisted FFT bump
    input.setup(gwstate);

    success = false; bool failed = false;
    try { run->run(gwstate, file_checkpoint, file_recoverypoint, logging);
          success = run->success(); file_progress.clear(); }
    catch (const TaskAbortException&) { if (!options.information_only) failed = true; }
    gwstate.done();

    if (success)       { primes++; composites = 0; if (!k_value.empty()) k_prime_found[k_value] = true; }
    else if (!failed)  { composites++; }
    if (failed && stop_error) Task::abort();
    if (Task::abort_flag()) break;
}
```

The two-`Logging` split is the key structural choice: `logging_batch` owns the *batch* progress file (`<batch>.param`, with `cur`/`primes`/`composites`) and prints the `"N of M: …"` lines, while a fresh `logging` per candidate owns that candidate's own `prst_<fingerprint>.param`/`.ckpt`/`.rcpt` and result lines — exactly as if it had been run standalone. That's what makes a batch resumable at two granularities: mid-candidate (the candidate's own checkpoint) and between-candidates (the batch `cur`).

Each candidate gets a **fresh `GWState` configured from the shared `Options`** — `options.configure(gwstate)` (`:332-333`) — so thread count / instruction set / safety margin parsed once at the top are shared, while FFT selection (`input.setup`) is per-candidate and can't leak into the next one. FFT-error state is per-candidate too: when a run hits repeated round-off errors, the `Task` error path bumps `gwstate.next_fft_count` and persists it as the `next_fft` param (`framework/task.cpp:140-141`); on resume, `:334-335` restores the bump from the candidate's own `.param` file. Note that `logging.file_progress(&file_progress)` is itself a **reset**: it clears the in-memory param map, re-reads it from the given file, and re-baselines the time accounting (`framework/logging.cpp:181-197`) — so the FFT bump and `time_total` a candidate sees come only from *its own* `.param` file, never from the previous candidate's run.

## 3. Field & method reference

The parser surface (`abc_parser.h`) is tiny; the driver's option set is the larger surface.

| `CandidateSource` member | Note |
|---|---|
| `size()` | candidate count; `0` is treated as "empty/unparseable" by the driver. |
| `get(index, out)` | fills `out`; `false` on out-of-range (the loop treats that as end-of-batch). Const + thread-safe by contract (header `:61`), though the driver calls it single-threaded. |
| `is_abc()` | true for ABC/ABCD/ABC2/NewPGen (carry k-values), false for raw. **Not called by the driver** (`batch.cpp` never references `is_abc()`). Whether `-stop on kprime` can act on a candidate is driven by that candidate's `k_value` being non-empty (`!k_value.empty()`, `batch.cpp:247`; the per-k record at `:361-363`), not by `is_abc()`. The two only *correlate* at source construction — raw sources are built with empty k-values and `_is_abc=false` (`abc_parser.cpp:965-968`), the header formats with k-values and `is_abc()==true` — but they're independent fields and can diverge (e.g. an ABC template with a leading `*` yields an empty `k_value` via `determine_k_variable`, yet `is_abc()` stays true). |
| `parse_batch_file(filename, logging)` | the factory (`abc_parser.cpp:879-969`): read lines → `detect_format` → dispatch → `CandidateSource` (or `nullptr` on failure). Raw is the fallback when no `ABC*`/NewPGen header matches. |

`batch_main` options (the `Config` DSL, `batch.cpp:48-138`) — a deliberate subset of `main()`'s, plus batch-only knobs:

| Option | Effect |
|---|---|
| `<file>` / `stdin` | the batch source (default-arg). `stdin` reads expressions line-by-line, unbounded. |
| `-stop on {error \| prime \| composites <n> \| primek}` | the stop conditions (§4); `kprime` is accepted as an alias of `primek` (`:121-122`). |
| `-newpgen {kn \| nk}` | NewPGen data-line column order: standard k-first (`kn`, default) or reversed n-first (`nk`) for merge scripts that emit `n k` (`:125`). |
| `-log [level] [batch <level>] [file <f>]` | separate verbosity for the batch log vs. per-candidate log. |
| `-info` | `print_info` each candidate instead of testing (continues unless `-fft info`). |
| `-trial` | trial-divide every candidate first; a found factor ⇒ "not prime", skip the full test. |
| `-ini <file>` | read more options from an ini file (`Config::parse_ini`). |
| `-t`, `-spin`, `-cpu`, `-fft`, `-check`, `-fermat`, `-order`, `-divides`, `-factors`, `-time`, `-d` | same meanings as in `main()`; parsed once into the shared `Options`, applied to each candidate's `GWState` via `options.configure` (`:333`). |

Note `Task::PROGRESS_TIME = 60` is set at entry (`:39`) — batch runs report progress less chattily than the 1-candidate default.

## 4. Lifecycle: file → source → loop → resume/stop

1. **`parse_batch_file`** (`abc_parser.cpp:879-969`): read all lines via `File::get_textreader`, `detect_format(line[0])`, then:
   - `FORMAT_ABCD` → `parse_abcd_file` → `VectorCandidateSource`.
   - `FORMAT_ABC2` → `parse_abc2_source` → `ABC2CandidateSource` (lazy).
   - `FORMAT_ABC` → parse the header as an `ABCTemplate`, expand each data line → `VectorCandidateSource`.
   - `FORMAT_NEWPGEN` → `parse_newpgen_file` (`:747-873`) → `VectorCandidateSource`.
   - else → **raw**: every line is an expression, no k-values, `is_abc() = false`.
2. **Resume bootstrap**: `batch_progress` (`<batch><suffix>.param`) is opened and `cur`/`primes`/`composites` are read back. `filename_suffix` is derived from `-order`'s fingerprint, `-divides` (`.div`), or `-fermat a` (`:188-193`) so two different bases over the same file don't collide.
3. **The loop** (§2): each iteration reports `cur`, saves batch progress, checks stop conditions, fetches the candidate, fast-paths small/trial numbers, then runs the full `Run`.
4. **Termination** (`:375-387`): if aborted **and** not stdin → save batch progress and return `PRST_EXIT_FAILURE` (so a re-run resumes); otherwise clear `batch_progress` and log `"Batch of N, primes: P, time: …"`.

**Stop conditions** (`-stop on …`):
- **`error`** — abort on a parse error *or* a task failure (`TaskAbortException` that wasn't info-only).
- **`prime`** — abort once any candidate is a (probable) prime. Checked at the *top* of the next iteration via the persisted `success` flag, so the prime's result is fully written first.
- **`composites <n>`** — abort after `n` **consecutive** composites; the counter resets to `0` on every prime (`:359-360`).
- **`primek`** (alias `kprime`) — once a prime is found for a given `k_value`, skip all later candidates sharing that `k` (`k_prime_found` map, `:247-252`). No-op for raw sources (no k-values).

**stdin mode** is the special case: no `source`, `total = 0`, the loop runs until an empty line or abort. It brackets the blocking `getline` with `time_total()`/`time_init` (`:231-233`) so the wall-clock spent waiting for input isn't billed to the candidate's timing.

## 5. The five input formats

All formats share line handling: `//` comments are stripped, blank/whitespace-only lines skipped (`is_skip_line`, `strip_comment`). `$a`…`$d` are the only variables.

**Raw** — no header recognized. Each line is passed verbatim to `InputNum::parse`. This is the fallback and the natural format for arbitrary, independent expressions (`100!-1`, `Phi(3,10^20)`, …); the ABC family instead applies one fixed template with numeric values substituted in. (A degenerate `ABC $a` template can also pass whole tokens through, so raw isn't strictly unique — but only single-token expressions survive, since `expand` splits the data line on whitespace.)

**ABC** — `ABC <template>` header, then one whitespace-delimited data row per candidate:
```
ABC $a*2^$b+1
5 100
7 100
```
→ `5*2^100+1`, `7*2^100+1`. The k-value is `$a` (the token before `*`).

**The sign convention (`+$c` vs a bare trailing `$c`).** Sieve tools write the tail of the template in two styles: the sign as a **template literal** (`ABC $a*2^$b+$c`, data value `1`) or a **bare `$c`** (`ABC $a*2^$b$c`) whose data value carries its own sign (`+1` / `-1`). Both work, because `expand` concatenates value tokens verbatim — signs included (`abc_parser.cpp:238-244`) — and `InputNum`'s tokenizer treats a leading `+`/`-` as an ordinary operator token (`framework/inputnum.cpp:395-400`): both styles produce the identical string `5*2^100+1`. The two conventions must not be mixed, though; there is no reconciliation code. Template `+$c` with a signed value produces `…+-1`/`…++1`, rejected by `InputNum::parse` as `"unexpected operator"` (`inputnum.cpp:889-897`); bare `$c` with an *unsigned* value glues the digits onto the exponent (`$b=100`, `$c=1` → `5*2^1001` — a syntactically valid but **wrong** number), caught only because the batch parses with `c_required=true` → `"C expected"` (`inputnum.cpp:959-965`). The k-variable extraction is unaffected either way — `determine_k_variable` only looks before the first `*` (`abc_parser.cpp:263-285`).

**ABCD** — delta-encoded, for runs of candidates that share most variables:
```
ABCD $a*2^$b+1 [5 100]
2
4
```
The `[…]` sets the initial accumulator per variable **and emits the first candidate** (`5*2^100+1`). Each subsequent line carries up to `num_vars` **deltas** added to the accumulators, then emits a candidate (`parse_abcd_file`, `:541-657`). So the rows above yield `5*2^100+1`, then `+2` on `$a` → `7*2^100+1`, then `+4` → `11*2^100+1`. Multiple `ABCD` header blocks may appear in one file, each resetting the accumulators. Accumulator arithmetic is `int64` with explicit overflow guards (`:618-627`).

**ABC2** — a template plus per-variable *range* definitions, expanded lazily as a Cartesian product:
```
ABC2 $a*2^$b+1
a: from 5 to 99 step 2
b: primes from 1000 to 2000
```
Each `<var>: …` line (`parse_abc2_var_line`, `:298-428`) is one of:
- **explicit list** `in { 5 7 11 }` (capped at `MAX_ABC2_IN_LIST_VALUES`, truncated if longer);
- **arithmetic range** `from X to Y [step S]` (default step 1) or `from X downto Y [step -S]`;
- **prime sieve** `primes from X to Y` (segmented sieve, span capped at `MAX_SIEVE_RANGE`; can't combine with `step`).

The candidate count is `∏ |values|`, checked against `MAX_ABC2_CANDIDATES` with a multiply-overflow guard (`:703-718`); `ABC2CandidateSource::get` then maps a flat index to one value per variable via precomputed strides.

**NewPGen** — the fixed-form sieve output of NewPGen (and the input format LLR reads). A colon-separated header, then `k n` pairs:
```
1000000000000:M:1:2:258
5 100
7 102
```
The header is `<sievelimit>:<char>:<chainlen>:<base>[:<mask>]` (`parse_newpgen_header`, `:74-84`); the form comes from the mask bits (`0x01` → `k*b^n+1`, `0x02` → `k*b^n-1`, per LLR's table, non-form flags `0x100`/`0x400` stripped first) or, for a 4-field header, from the char code (`P`/`M`). The rows above yield `5*2^100-1`, `7*2^102-1`; the k-value is the `k` column. Everything NewPGen can express beyond the two single-test forms — twins, Sophie-Germain, CC chains, primorials, `+5`/`+7`, AP, dual — is **warned about and the whole block skipped** (`parse_newpgen_file`, `:747-873`). Multiple header blocks per file are fine (each resets the form/base). `-newpgen nk` flips the data columns for merge scripts that emit `n k`.

## 6. Pitfalls

- **The small-number and `-trial` fast paths bypass the prime/composite accounting.** A candidate with `bitlen() ≤ 40` (`:280`) or a `-trial` factor hit (`:301-313`) prints its result and `continue`s *before* the `primes`/`composites`/`k_prime_found` bookkeeping (`:356-369`). So those candidates don't increment the composite streak, don't count toward `-stop on composites`, and a small prime won't arm `-stop on primek`. Surprising if a batch is all small numbers.
- **`-stop on prime` stops one iteration late, by design.** `success` is checked at the loop top (`:213`), not right after the run, so the prime's result line and checkpoint cleanup complete first. Don't read the abort as "stopped mid-prime."
- **The abort-based stop conditions exit with `PRST_EXIT_FAILURE` (1), even on success.** `-stop on prime`, `on composites`, and `on error` all call `Task::abort()`; the post-loop check (`:376`) maps a set abort flag (when not stdin) to `progress_save()` + `return PRST_EXIT_FAILURE`. So a *successful* stop-on-prime returns exit 1, not 0 — a script keying off the exit code reads it as failure. Only a batch that runs to natural completion returns `PRST_EXIT_NORMAL` (0). (`-stop on primek` is the exception: it `continue`s to skip candidates and never aborts, so it doesn't force the failure exit.)
- **`-stop on composites` counts *consecutive* composites, not total.** Any prime resets the counter to 0 (`:359-360`). A batch that alternates prime/composite never trips it.
- **ABC data values and template signs must agree.** `+$c` templates want unsigned values, bare-`$c` templates want signed ones; a mismatch either hard-fails (`…+-1` → "unexpected operator") or glues digits onto the exponent (`…^1001` → "C expected") — see §5. The rejection happens per data line, so under the default (no `-stop on error`) a mismatched file quietly skips *every* candidate.
- **`detect_format` requires a trailing space for the ABC family.** `"ABC "`, `"ABCD "`, `"ABC2 "` (`:86-93`). A header line `ABC2\t…` (tab, no space) or `ABCD[…]` glued to the bracket falls through — usually to **raw** (a NewPGen misread is unlikely: the line would have to scan as ≥ 4 colon-separated numeric fields, `:95-103`), where each subsequent line is then mis-parsed as a standalone expression.
- **Raw is the silent fallback.** `parse_batch_file` never fails on "unknown format" — anything without a recognized `ABC*`/NewPGen header becomes a raw source. A typo'd `ABBC` header means the header line itself becomes candidate #1 (and fails `InputNum::parse`).
- **NewPGen supports only `k*b^n±1`.** Twin/SG/CC-chain/primorial/AP/dual masks make the whole block a warning + skip (`:796-819`); a file that is *all* unsupported blocks parses to zero candidates and the batch reports "no supported candidates."
- **`filename_suffix` only disambiguates `-order`/`-divides`/`-fermat a`.** Two batch runs of the same file with *different* `-check`/`-factors` options share the same `.param`/`.ckpt` files. Resuming one over the other can read a stale `cur`. Use distinct working directories if running variants concurrently.
- **ABC2 limits silently truncate or refuse.** An `in { … }` list over 1000 entries is truncated (warning only); a Cartesian product or prime span over the caps makes `parse_abc2_source` return `nullptr` → the whole batch is "empty or could not be parsed." Check the warnings.

## 7. Quick reference

| Format | Header (line 1) | Body | Source impl | k-value |
|---|---|---|---|---|
| raw | *(none recognized)* | one expression per line | `VectorCandidateSource` | none |
| ABC | `ABC $a*2^$b+1` | whitespace-delimited value rows | `VectorCandidateSource` | token before `*` |
| ABCD | `ABCD $a*2^$b+1 [init…]` | per-variable deltas; multiple blocks ok | `VectorCandidateSource` | token before `*` |
| ABC2 | `ABC2 $a*2^$b+1` | `<var>: in{…}` / `from…to…[step]` / `primes from…to` | `ABC2CandidateSource` (lazy) | token before `*` |
| NewPGen | `<limit>:<P\|M>:<len>:<base>[:<mask>]` | `k n` rows (`n k` with `-newpgen nk`); multiple blocks ok | `VectorCandidateSource` | the `k` column |

| You want to… | Where |
|---|---|
| Add a batch option | the `Config` chain in `batch_main` (`batch.cpp:48-138`) |
| Change format detection | `detect_format` (`abc_parser.cpp:86-106`) |
| Add a new ABC2 range keyword | `parse_abc2_var_line` (`abc_parser.cpp:298-428`) |
| Support another NewPGen form (twin, SG, …) | the mask/char dispatch in `parse_newpgen_file` (`abc_parser.cpp:790-825`) |
| Resume a batch | re-run the same command; `cur`/`primes`/`composites` come from `<batch><suffix>.param` |
| Stream candidates from another process | pipe into `PRST -batch stdin …` |
| Test arbitrary expressions (not a sieve template) | use a raw file (one expression per line) |

## 8. Open questions / non-coverage

- **`Config` DSL internals.** The option grammar (`group`/`exclusive`/`value_code`/`parse_ini`) is `framework/config.{h,cpp}`, shared with `main()`; this doc only lists `batch_main`'s resulting options. The DSL itself is now documented in **`config-dsl.md`**.
- **Thread-safety of `CandidateSource`.** The header promises const-`get` thread-safety (`abc_parser.h:57`) but the driver is single-threaded. If batch-level parallelism is ever added, `ABC2CandidateSource::get` (which builds throwaway vectors per call) is already safe, but the `Logging`/`Progress` bookkeeping is not — see `logging-and-progress.md` §12.
- **The exact ABCD multi-block / per-variable-advance semantics.** `parse_abcd_file` advances accumulators one token-set per data line and supports multiple `ABCD` headers; the precise interaction with sieve tools' emitted deltas is reproduced faithfully from the code but not cross-checked against srsieve2's writer. Verify against a real sieve file if a delta-decoding bug is reported.
- **`File::get_textreader` / `read_buffer`.** The line-reading path is the `File` abstraction's text mode; covered (as on-disk I/O) by the framework's `state-serialization.md`. Here it's just "read the lines."
