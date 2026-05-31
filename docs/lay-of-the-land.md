# PRST — Lay of the Land

A high-level orientation to the codebase for contributors new to PRST. For deeper dives see the per-subsystem docs in this folder.

## What PRST is

PRST is a C++ command-line utility for **primality testing of very large numbers**. Written by Pavel Atnashev (`patnashev/prst` on GitHub), built on George Woltman's GWnum multiplication library. It is best at systematic searches of numbers in popular forms — Proth, Thabit, generalized Fermat, factorials, primorials, arbitrary candidates — and is typically driven by a sieving stage upstream of it. (Mersenne numbers go to GIMPS instead.)

It implements:
- **Fermat** probabilistic test (any number, not definitive)
- **Proth** test of `k*2^n+1`
- **Pocklington** test of `k*b^n+1`, `n!+1`, `n#+1`
- **Morrison** test of `k*b^n-1`, `n!-1`, `n#-1`
- **LLR** as a special case of Morrison

Plus a verification-grade proof system (`-proof save/build/cert`) that lets a tester emit a certificate which any third party can verify cheaply, an order-finding mode (`-order`), and a batch driver (`-batch`).

Two optional build flavors are gated by `#define`s in `src/version.h`:
- `BOINC` — distributed-computing wrapper
- `NETPRST` — networked client

When defined, each adds a runtime subcommand to `main()`'s option parser (`-boinc` at `prst.cpp:199-201`, `-net` at `prst.cpp:202-204`) which `exit()`s into `boinc_main` / `net_main` (in `boinc.cpp` / `net.cpp`). They reuse the same test machinery; they don't replace `main()`.

## Repo layout

```
prst/                                ← repository root (fork of patnashev/prst)
├── README.md
├── sample.ini                       ← example config
├── docs/                            ← this folder
├── src/                             ← PRST-specific code (~10k LoC)
│   ├── prst.cpp / prst.h        ← entry point, option parsing, Run::create dispatcher
│   ├── fermat.cpp / fermat.h    ← Fermat / Proth (also base for Pocklington)
│   ├── pocklington.cpp / .h     ← Pocklington + PocklingtonGeneric
│   ├── morrison.cpp / .h        ← Morrison + MorrisonGeneric
│   ├── order.cpp / .h           ← multiplicative order
│   ├── proof.cpp / .h           ← proof certificates (save/build/verify)
│   ├── exp.cpp / .h             ← exponentiation tasks (CarefulExp, MultipointExp …)
│   ├── lucasmul.cpp / .h        ← Lucas-chain tasks (LucasVMul, LucasUVMul …)
│   ├── batch.cpp                ← -batch mode, parses ABC files
│   ├── abc_parser.cpp / .h      ← ABC/ABCD/ABC2 batch-file parser (newer split)
│   ├── boinc.cpp / boinc.h      ← BOINC entrypoint
│   ├── net.cpp / net.h          ← NETPRST client
│   ├── testing.cpp / testing.h  ← built-in self-tests (`-test`)
│   ├── support.cpp / .h         ← LLR2File — LLR2 checkpoint-format compatibility (File subclass)
│   ├── version.h                ← version + build feature flags
│   ├── win64/  linux64/  mac64/ ← per-platform build files
│   └── test.data                ← canned inputs for self-tests
└── framework/                   ← submodule → patnashev/arithmetic
    ├── arithmetic/              ← Giant, GWArithmetic, lucas, edwards, montgomery, poly …
    ├── gwnum/                   ← prebuilt GWnum library (Woltman's FFT-based bignum)
    ├── gmp/                     ← prebuilt GMP for Windows
    ├── bow/                     ← Big Object Writer (state serialization)
    ├── inputnum.cpp / .h        ← number parser + form classification
    ├── logging.cpp / .h         ← Logging / SubLogging / Progress
    ├── task.cpp / .h            ← Task / InputTask base classes
    ├── file.cpp / .h            ← File abstraction with checkpoints
    ├── config.cpp / .h          ← option-parser DSL used in prst.cpp
    ├── container.cpp / .h       ← FileContainer (proof packs)
    └── md5.c / md5.h            ← hashing
```

Mental model: **`src/` is "what test do we run and how do we orchestrate it?"; `framework/` is "the shared bignum + arithmetic + checkpointing primitives every test reuses."** The framework lives in a separate repo (`patnashev/arithmetic`) and is consumed by both PRST and other Atnashev-family tools.

## End-to-end execution flow

A typical invocation `prst "30006!4-1" -d` flows:

1. **`main()` in `src/prst.cpp:43`** — install signal handlers, set up `GWState`, declare `Options`.
2. **Option parsing** — uses the `Config` DSL (`framework/config.h`) to parse argv into `Options`, `gwstate`, `proof_op`, log level, etc. (`prst.cpp:82` onward.)
3. **Input parsing** — `InputNum::parse` (`framework/inputnum.cpp`) classifies the candidate as one of `KBNC` / `FACTORIAL` / `PRIMORIAL` / `GENERIC` and stores k, b, n, c, factors.
4. **Trial-division shortcut** — small inputs (`bitlen <= 40`) or `-trial` are handled inline at `prst.cpp:292`.
5. **Progress file wiring** — a `File` for `prst_<fingerprint>*.param` is opened so progress survives restarts (`logging.file_progress(...)` at `prst.cpp:334`).
6. **`Run::create()` (`prst.cpp:444`)** — the dispatcher, in code order:
   - `-order` → `Order`
   - `-proof cert` → reuse the `Proof` instance directly
   - forced `-fermat`, generic numbers, or `|c| ≠ 1` → `Fermat::FERMAT`
   - Proth-form `k*2^n+1` with `k < 2^n` → `Fermat::PROTH`
   - call `expand_factors()`; if not half-factored → `Fermat::AUTO` with a "Not enough factors" warning
   - `c == 1` and half-factored → `Pocklington` (when `n > 10` and `factors.size() < 10`) else `PocklingtonGeneric`
   - `c == -1` and half-factored → `Morrison` (when `n > 10` and `factors.size() < 10`) else `MorrisonGeneric`
   - otherwise → `nullptr` (failure)
7. **Optional proof wrapping** — if `-proof save|build` is requested, the `Run` is wrapped in a `Proof` that delegates the inner Fermat test (`prst.cpp:358-381`).
8. **`run->run(gwstate, file_checkpoint, file_recoverypoint, logging)`** — the polymorphic test entrypoint at `prst.cpp:419`. The selected test class drives one or more `Task` instances, writes checkpoints to `prst_<fingerprint>.ckpt`/`.rcpt`, and emits result lines via `logging.result(...)`.
9. **Exit code** — `PRST_EXIT_PRIMEFOUND` (2) / `PRST_EXIT_NORMAL` (0) / `PRST_EXIT_FAILURE` (1), defined in `src/prst.h`.

## The Run hierarchy

Introduced in upstream commit `4764aa9 Refactoring, base class for all tests.` All tests now derive from a common `Run` (`src/prst.h:52`):

```
Run                         (abstract base — name, fingerprint, success/prime flags, run())
├── Fermat                  (Fermat / Proth / Pocklington-base / AUTO; uses MultipointExp)
│   └── Pocklington         (small-factor Pocklington; reuses Fermat::run)
├── PocklingtonGeneric      (FactorTree-driven; runs sub-tasks under a SubLogging)
├── Morrison                (small-factor Morrison; uses LucasVMul)
├── MorrisonGeneric         (FactorTree-driven; runs sub-tasks under a SubLogging)
├── Order                   (multiplicative order; MultipointExp + CarefulExp)
└── Proof                   (orchestrates a wrapped Fermat for -proof save/build/cert)
```

Two patterns repeat:
- **Small-factor variants** (`Morrison`, `Pocklington`) work directly on the parent `Logging`. They are picked when `n > 10` and `factors.size() < 10` (`prst.cpp:480, 490`).
- **`*Generic` variants** (`MorrisonGeneric`, `PocklingtonGeneric`) build a `FactorTree`, walk it, and run sub-tasks under their own `SubLogging _logging`. This is where sub-task time accounting happens, and it is a subtle area (the inner `SubLogging` progress must be propagated to the parent).

`Run::create` (a static factory in `prst.cpp:444`) is the single source of truth for which class handles which form.

## The Task layer

Tests don't do the heavy math themselves — they construct one or more `Task` subclasses from the framework + `src/exp.h` / `src/lucasmul.h`, configure them, and call `task->run()`. The Task base lives in `framework/task.h:61`:

- `Task` — abstract; manages `_state`, checkpointing cadence (`MULS_PER_STATE_UPDATE`, `DISK_WRITE_TIME`, `PROGRESS_TIME`), the global `_abort_flag` (signal handlers call `Task::abort()` from `prst.cpp:38`), and the `commit_setup`/`commit_execute` lifecycle.
- `InputTask : Task` — adds `_input`, `_timer`, error-check toggles. This is the more common base.

Concrete tasks of interest:
- `CarefulExp`, `MultipointExp`, `BaseExp` — exponentiation variants with strong-error-check (Gerbicz/Gerbicz-Li) options. Defined in `src/exp.h`.
- `LucasVMul`, `LucasVMulFast`, `LucasUVMul` — Lucas-sequence multipliers. Defined in `src/lucasmul.h`.
- `Product` — accumulates a product of giants. Used in MorrisonGeneric/PocklingtonGeneric.

A test class typically:
1. Constructs its tasks in its constructor (using `add_stage` to register expected work for progress reporting).
2. Resumes from `file_recoverypoint` if a state file exists.
3. Loops: call `task->run()`, inspect `task->result()`, advance the algorithm.
4. Writes the result line via `logging.result(...)` and `logging.result_save(...)`.

## Cross-cutting infrastructure

### `Logging` / `SubLogging` / `Progress` — the hot zone

(Detailed write-up: see `logging-and-progress.md` in this folder.)

- `Logging` (`framework/logging.h:50`) owns one `Progress`, the result file, the factor file, the log file, and the progress checkpoint file. All `info/warning/error/result` flow through `report()`.
- `SubLogging` (`framework/logging.h:108`) wraps a parent. Used when a test runs an inner sub-task that should see its own progress, but whose results, params, factors, and progress saves still bubble up to the parent. **Its constructor calls `progress().set_parent(&parent.progress())`** — an easy detail to overlook.
- `Progress` tracks staged costs (`add_stage`), current stage progress, and time. It maintains:
  - `_time_total` — durable; only `update()` and `time_init()` change it.
  - `_time_stage` — transient; reset to 0 by every `update()` call on this Progress; receives `+= elapsed` walks from descendants.
  - `time_total() = _time_total + _time_stage`.
  - `update()` walks up `_parent` and pokes `_time_stage` / `_cur_progress` on each ancestor.

Non-obvious invariants worth pinning:
- `_time_stage` **is meant to be transient** between successive `update()` calls. Anyone relying on it persisting after a sibling `update()` is wrong.
- Time only reaches `_time_total` of an outer Progress if `update()` is called on that outer Progress (capturing the wall-clock since its previous `update`). Sub-progress updates do NOT increment ancestors' `_time_total`.

### `InputNum` — the parsed candidate

`framework/inputnum.h:10`. Stores k/b/n/c plus factor lists and pretty-printed text. Exposes `type()` (`KBNC`/`FACTORIAL`/`PRIMORIAL`/`GENERIC`), `value()`, `bitlen()`, `fingerprint()` (32-bit hash used for filename uniqueness), `factors()`, `cofactor()`, `is_half_factored()`. Most of the dispatcher's branching keys off this.

### `File` — checkpoint storage

`framework/file.h`. Filename + 32-bit fingerprint header + appid (`File::FILE_APPID = 4` for PRST). Has children for sub-state. The progress file (`*.param`) and checkpoint files (`*.ckpt`, `*.rcpt`) all use this. `LLR2File` is a variant with a different on-disk format used for some compatibility cases (`prst.cpp:336`).

### `GWState` — the math runtime

Owned by `main()`. Configures GWnum: thread count, spin threads, instruction set (SSE2/AVX/FMA3/AVX512F), FFT size hints, safety margin, fingerprint. `input.setup(gwstate)` picks the actual FFT and reports it via `gwstate.fft_description`.

### `Options` — user-facing knobs

`src/prst.h:10`. Bag of `std::optional`s populated by the `Config` DSL. Test classes consult fields like `Check`, `CheckStrong`, `StrongCount`, `FermatBase`, `AllFactors`, `ProofPointFilename`. `Run::create` may also set `options.maxmulbyconst` as a side-effect.

## Build

- Windows: `Build commands.txt` invokes MSBuild on `src/win64/PRST.sln` (Release/x64, Win SDK 10.0.26100.0). `BOINC`/`NETPRST` builds are separate solution configs.
- Linux: `src/linux64/Makefile` (and `Makefile.boinc`).
- macOS: `src/mac64/Makefile`.

GWnum + GMP are prebuilt static libs shipped in the framework submodule (`framework/gwnum/{linux64,mac64,win64}`, `framework/gmp/win64`). Bumping GWnum means bumping the framework submodule pointer.

## Working with this codebase — practical notes

1. **Read the framework before touching it.** The framework is a shared library used by other tools; its surface area is tighter than `src/`. Adding a method to `framework/logging.h` is a much higher-cost change than modifying a caller in `src/morrison.cpp`.
2. **`Run::create` is the dispatcher.** Any new test mode plugs in there. Any change to which class handles which input form lives there.
3. **Result lines are user-visible contract.** The `logging.result(...)` strings (`is prime!`, `is not prime`, `time: %.1f s.`) are what mersenneforum users / BOINC parsers consume — don't change wording lightly.
4. **Checkpoint format is on-disk contract.** `TaskState::read/write` and the `File` fingerprint scheme define resumability. Bumping any of them needs a version bump and ideally backward-compat reads.
5. **Single-repo fixes preferred.** Coordinated `prst` + `framework` PRs are reviewable but unloved by upstream. Default to "can I express this with the existing framework API?" before reaching for a new one.
6. **Bug reports come from mersenneforum.org.** Search the forum threads for the exact result-line snippet a user is complaining about — that's typically the fastest path to the failing code path.

## Companion deep-dives in this folder

- **`logging-and-progress.md`** — `Logging` / `SubLogging` / `Progress`, the parent walk in `update()`, persistence, common patterns, pitfalls.
- **`task-lifecycle.md`** — `TaskState` / `Task` / `InputTask`, the restart loop in `run()`, the cadence in `on_state()`, error-correction handshake, concrete subclass tour, pitfalls.
- **`proof-system.md`** — the `-proof save / build / cert` pipeline; `Proof` as both `Run` and Fermat-wrapper; `ProofSave` / `ProofBuild` tasks; the binary product-tree schedule; roots-of-unity and security-seed defenses; the `.pack` container.
- **`run-hierarchy.md`** — Fermat / Pocklington(Generic) / Morrison(Generic) / Order at the call-graph level; the `Run::create` dispatcher; the small-factor vs. `*Generic` boundary; `FactorTree` construction and walk; the result-line contract.
- **`inputnum-parsing.md`** — `InputNum::parse` classification (`KBNC`/`FACTORIAL`/`PRIMORIAL`/`GENERIC`); the `(<num>)/F` cofactor-divisor form; `Phi`/`Quad`/`Hex` and auto-detected algebraic factoring; the negative-`Giant` factorial/primorial encoding; `mod`/fingerprinting; `is_half_factored`.
- **`abc-batch.md`** — `-batch` as an alternate entry point; `CandidateSource` + `parse_batch_file`; the raw / ABC / ABCD / ABC2 formats and `ABCTemplate`; the `batch_main` driver loop; two-level progress/resume (`cur`/`primes`/`composites`); the `-stop on {error,prime,composites,kprime}` conditions.
- **`boinc-and-net.md`** — the `BOINC` / `NETPRST` build flavors; `boinc_main` / `net_main` as `exit()`-into-sibling entry points; `BoincLogging` (control hooks: `state_save_flag`/`heartbeat`/trickles) vs. `NetLogging` (data hooks); `NetFile`/`LLR2NetFile` checkpoint shipping over HTTP; the reduced inline dispatcher; **`net.cpp` is stale vs. the run() refactor**.
- **`state-serialization.md`** — the `Writer`/`Reader`/`File` on-disk format (magic + appid + `TYPE` + version + fingerprint header); the `TYPE` registry `{1,2,3,4,6,8,9,10,11}` and their `TaskState` owners; `Giant`/`SerializedGWNum` encoding; atomic write + `.md5` sidecar; `unique_fingerprint`; the `LLR2File` compatibility munging (byte 4 = appid).
- **`arithmetic-foundation.md`** — the bignum substrate: the `Giant` heap model (`giant_struct`, ABI-compatible with gwnum's C `giant`) and its `Giants`/`GMP`/`GW` backends; `GWState::setup` (the three FFT-config paths, fingerprint, `known_factors` reduction); `GWArithmetic`/`GWNum`; `SerializedGWNum` (the checkpoint bridge); the `ReliableGWArithmetic` round-off escalation (`restart_flag`/`failure_flag`/`suspect_ops`).
- **`config-dsl.md`** — the `Config` option-parser DSL: the runtime `ConfigObject` tree vs. the CRTP `*Setup` builder; the `delim` matching convention (glued / next-token / literal); the builder verbs (`check`/`value_*`/`group`/`exclusive`/`list`/`on_check`/`default_code`); the `parse_args` walk + group fixpoint + exclusive first-match; `.ini` parsing.
- **`test-harness.md`** — `PRST -test`: the `Test`/`DeterministicTest` classes and the `test.data` vector arrays (`res64`/`cert64` oracle); `Test::run`'s SAVE→BUILD→CERT cross-check (incl. the cross-FFT consistency check); the subsets; `RootsTest` (proof-forgery defenses) and `ABCParserTest`; the separate standalone `arithmetic/test.cpp`.
- **`exponentiation-algorithms.md`** — the exp/Lucas `InputTask`s that do the modular work: `BaseExp`/`MultipointExp` (point schedule, sliding window, smooth-vs-non-smooth) and `StrongCheckMultipointExp`/`LucasUVMul` (the Gerbicz / Gerbicz-Li block check — `L≈√iters`, the `D` check-accumulator, recovery-point rollback); how `Fermat::run` picks the task; the `State` TYPEs (1/2/8, 9/10/11).
- **`curves-and-polynomials.md`** — the specialized arithmetic *beside* the bignum layer: the full (`GroupArithmetic`, NAF) vs. differential (`DifferentialGroupArithmetic`, DAC ladder) group split; the five concrete arithmetics; `get_NAF_W`/`precomputed_DAC_S_d` chain construction. Key liveness fact: only `LucasV`/`LucasUV` (Morrison) is a **live** PRST consumer — the Edwards/Montgomery ECM curves are dormant (commented-out caller) and `PolyMult` has no live caller (the proof tree multiplies giants, not polynomials). An interface+usage map (the ~42k-line `group.cpp` formulas are fenced to the math doc).
- **`math-and-theorems.md`** — the *why*: Fermat PRP vs. the deterministic Proth (iff) and Pocklington/BLS (`N−1`, `F>√N`) and Morrison/BLS-Theorem-14 (`N+1`, Lucas) criteria; the half-factored threshold; Gerbicz / Gerbicz-Li error checking; IBDWT multiplication; the theorem→code map. A reading guide tying each test's correctness to its code and to the authoritative references (BLS 1975, eprint 2023/195, the `mult_*.pdf`s) — proofs cited, not reproduced.
- **`container-format.md`** — the `.pack` bundle format (live via `-proof … pack`): the newline-delimited-JSON-record on-disk layout (magic header; `[size, stream_id, codec]` chunks; `stream_id 0` index chunks of `{file,stream,size,offset,md5}`; the append/override model); `Packer`/`FileContainer`/`FilePacked`; codecs, MD5, and the `container_error` recovery codes. Grounded in the `FileContainer::open` parser; the writer byte-framing is fenced.

## Coverage status

**Every *in-scope* subsystem has a deep-dive** (the fifteen in "Companion deep-dives" above); the only things without one are explicitly scoped out (the third-party prebuilts, `bow.cpp`'s BOINC-only build glue, and the fenced formula/byte details — all listed below). A durable caveat worth keeping in mind: **a "complete" claim has to survive a grep one level deeper, and lower-traffic, lower-level code is *more* worth documenting, not less.**

What's *intentionally* left without its own doc (not gaps — scoped-out by design):

> **In scope vs. external.** `framework/arithmetic/` is patnashev's own source (the `patnashev/arithmetic` submodule) — documented in `arithmetic-foundation.md` + `curves-and-polynomials.md`. By contrast `framework/gwnum/` (Woltman's prebuilt FFT bignum) and `framework/gmp/` are **prebuilt third-party libraries**: PRST depends on their *interface* (`gwnum.h`, GMP) but owns no source, so they stay black boxes by nature — understood at the API level (mostly via `arithmetic.cpp`'s wrappers). `framework/bow/bow.cpp`'s implementation is BOINC-only build glue, relevant only under the `BOINC` build, so it stays out of scope (its "Big Object Writer" label in the repo-layout section is itself unverified — see `boinc-and-net.md` §8).

- **`main()`'s top half** — option parsing, the `bitlen ≤ 40` trial-division shortcut, and progress-file wiring before `Run::create` — is walked step-by-step in "End-to-end execution flow" above; that counts as covered.
- **The full mathematical proofs** (BLS Theorem 14, the Gerbicz-Li soundness bound, the IBDWT error bound) live in the cited references; `math-and-theorems.md` is the reading guide to them, not a re-derivation.
- **`group.cpp`'s ~42k lines of curve/poly formulas** and `container.cpp`'s writer byte-framing/codec internals are fenced in their docs (`curves-and-polynomials.md`, `container-format.md`) — interface + on-disk layout documented, formula/byte details pointed to.

A caution for maintainers: the docs reflect the code at the time of writing and several embed verified line numbers and "this is the only live caller" claims — re-grep before trusting any such claim.

## Cross-cutting open questions (parked)

These came up while writing the existing docs and didn't fit into one subsystem:

- **Thread-safety boundary.** PRST's threads live inside GWnum's FFT layer; everything above is single-`Task` on the main thread. If multi-`Task` parallelism is ever added (e.g. parallel inner factors in `MorrisonGeneric`), the entire `Logging`/`Progress` subsystem needs locking. See `logging-and-progress.md` §12 and `task-lifecycle.md` Pitfall E.
- **`ReliableGWArithmetic` semantics.** ~~Treated as a black box in `task-lifecycle.md` §8.~~ **Resolved:** the round-off escalation (`restart_flag`/`failure_flag`/`suspect_ops`/`op`) is documented in `arithmetic-foundation.md` §5.
- **`_smooth` exponentiation path.** ~~`BaseExp::_smooth` toggles a different math path.~~ **Resolved:** the smooth (`b^n` by repeated squaring / windowed powering) vs. non-smooth (sliding-window over the full exponent) split is documented in `exponentiation-algorithms.md` §1, §2.

## How to write the next deep-dive

The template is `docs/logging-and-progress.md` (and now `docs/task-lifecycle.md`). Sections, in this order:

1. One-paragraph intro + source file list + cross-references.
2. Class hierarchy / data-model walkthrough.
3. Annotated source for the central function/method (verbatim from the codebase, with inline comments).
4. Field/method reference, every line annotated — no bare lines.
5. Lifecycle / sequence section.
6. Common patterns.
7. Pitfalls, anchored to real bugs where possible.
8. Quick-reference table.
9. Open questions / explicit non-coverage, with pointers to the relevant other doc(s).

When a new deep-dive lands, add its bullet to "Companion deep-dives in this folder" above, and remove anything it now covers from a sibling doc's "open questions."
