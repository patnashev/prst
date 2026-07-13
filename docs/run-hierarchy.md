# Run hierarchy — deep dive

Every primality test PRST can run is a subclass of `Run`. `main()` parses the candidate, `Run::create()` picks the right subclass for the candidate's *form*, and `run->run(...)` drives the test. This doc is the call-graph-level tour of that layer: what each subclass's `run()` actually does step by step, why some forms get a "small-factor" class and others a `*Generic` one, how `FactorTree` is built and walked, and the exact user-visible result lines each test emits.

The hierarchy was introduced by upstream commits `4764aa9 Refactoring, base class for all tests.` and `fe25a5a Refactoring of run().` — it is the freshest part of `src/` and the least documented. The dispatch logic in `Run::create` is denser than it looks; any new test mode or any change to which class handles which input form plugs in here.

Source files:
- `src/prst.h:75-100` (the `Run` base class), `prst.h:10-71` (`Options`)
- `src/prst.cpp:463-529` (the `Run::create` dispatcher)
- `src/fermat.{h,cpp}` (`Fermat`: FERMAT / PROTH / AUTO / POCKLINGTON-base)
- `src/pocklington.{h,cpp}` (`Pocklington`, `PocklingtonGeneric`, and `FactorTree`)
- `src/morrison.{h,cpp}` (`Morrison`, `MorrisonGeneric`)
- `src/order.{h,cpp}` (`Order` and `FermatDivisor`)

Prereqs: `task-lifecycle.md` (every `run()` here constructs `InputTask`s and calls `task->run()` — the restart/checkpoint loop lives there), `logging-and-progress.md` (the `*Generic` tests run their inner tasks under a `SubLogging`, and `add_stage`/`next_stage`/`update` are the progress contract), `proof-system.md` (the `Proof` wrapper that sits *above* `Fermat` for `-proof save|build`).

## 1. The hierarchy at a glance

```
Run                         (abstract base — name, fingerprint, success/prime flags, run())
├── Fermat                  (Fermat / Proth / Pocklington-base / AUTO; uses MultipointExp)
│   └── Pocklington         (small-factor Pocklington; reuses Fermat::run)
├── PocklingtonGeneric      (FactorTree-driven; runs sub-tasks under a SubLogging)
├── Morrison                (small-factor Morrison; uses LucasVMul)
├── MorrisonGeneric         (FactorTree-driven; runs sub-tasks under a SubLogging)
├── Order                   (multiplicative order; MultipointExp + CarefulExp)
│   └── FermatDivisor       (-divides: F/GF/xGF divisibility search; reuses Order's smooth tasks)
└── Proof                   (orchestrates a wrapped Fermat for -proof save/build/cert)
```

`Proof` is covered in `proof-system.md`; it is the one subclass that *contains* a `Run` rather than being a leaf. Everything else here is a leaf except `Order`, which `FermatDivisor` derives from.

Two structural facts that are easy to misread:

- **Two tests derive from another test.** `Pocklington : public Fermat` (`pocklington.h:6`) reuses Fermat's exponentiation machinery and `Fermat::run` for the initial probable-prime stage; `FermatDivisor : public Order` (`order.h:49`) reuses Order's smooth-task machinery via the protected `Order(const char*, InputNum&, Options&)` constructor and `create_smooth_task`. The rest (`PocklingtonGeneric`, `Morrison`, `MorrisonGeneric`, `Order`) derive directly from `Run` and build their own tasks from scratch.
- **`Fermat` declares two `run` overloads.** The `Run::run` override (`fermat.h:23`) just forwards to a second overload that takes an extra `Proof*` (`fermat.h:24`). `Pocklington` overrides *that* second overload (`pocklington.h:11`). This is the proof-callback protocol described in `proof-system.md` §intro — don't "simplify" it away.

Constructor signatures (note which take a `Proof*`):

| Class | Constructor | Takes `Proof*`? |
|---|---|---|
| `Fermat` | `Fermat(int type, InputNum&, Options&, Logging&, Proof*)` | yes |
| `Pocklington` | `Pocklington(InputNum&, Options&, Logging&, Proof*)` | yes (forwards to `Fermat`) |
| `PocklingtonGeneric` | `PocklingtonGeneric(InputNum&, Options&, Logging&)` | no |
| `Morrison` | `Morrison(InputNum&, Options&, Logging&)` | no |
| `MorrisonGeneric` | `MorrisonGeneric(InputNum&, Options&, Logging&)` | no |
| `Order` | `Order(InputNum&, Options&, Logging&)` | no |
| `FermatDivisor` | `FermatDivisor(InputNum&, Options&, Logging&)` | no |

Proofs are only supported on the Fermat-family path; the `*Generic` and `Morrison`/`Order` paths take no `Proof*` (and `PocklingtonGeneric::run` aborts if it ever needs to restart while a proof is active — see §4.5/§7).

## 2. `Run` — the base class

`prst.h:52-77`, verbatim:

```cpp
class Run
{
public:
    Run(InputNum& input_, Options& options) : input(input_), _options(options) { }
    Run(const char* name, InputNum& input_, Options& options) : _name(name), input(input_), _options(options) { }
    virtual ~Run() { }

    const std::string& name() { return _name; }
    uint32_t fingerprint() { return _fingerprint; }
    bool success() { return _success; }
    bool prime() { return _prime; }
    std::string& res64() { return _res64; }

    virtual void run(arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging) = 0;

    static Run* create(InputNum& input, Options& options, Logging& logging, Proof* proof = nullptr);

protected:
    InputNum& input;
    Options& _options;
    std::string _name;
    uint32_t _fingerprint = 0;
    bool _success = false;
    bool _prime = false;
    std::string _res64;
};
```

Field/method reference — every member:

| Member | Role |
|---|---|
| `input` | the parsed candidate (`InputNum&`, a reference — not owned). All dispatch and exponent construction reads `input.factors()`, `input.c()`, `input.type()`, `input.value()`, etc. |
| `_options` | user knobs (`Options&`, see below). Tests read `CheckStrong`, `StrongCount`, `FermatBase`, `AllFactors`. |
| `_name` | human-readable test name (`"Fermat test"`, `"Proth test"`, `"Morrison (LLR) test"`, …). Set in the ctor, sometimes refined mid-ctor once the exact form is known. |
| `_fingerprint` | 32-bit ID used to namespace checkpoint files. Usually `input.fingerprint()` or `File::unique_fingerprint(input.fingerprint(), <salt>)`. |
| `_success` | the candidate passed the probable-prime stage (Fermat/Lucas congruence held). Not yet a proof of primality. |
| `_prime` | the candidate is *proved* prime (Pocklington/Morrison factor conditions met, or Proth result was 0). |
| `_res64` | RES64 hex of the failing residue (or the GCD'd factor) — populated only on the not-prime / factor-found paths. |
| `name()` / `fingerprint()` / `success()` / `prime()` / `res64()` | accessors `main()` reads after `run()` returns (e.g. to pick the exit code). |
| `run(gwstate, file_checkpoint, file_recoverypoint, logging)` | the one pure-virtual entrypoint. `file_checkpoint`/`file_recoverypoint` are the `.ckpt`/`.rcpt` `File`s; tests add per-base child files under them. |
| `create(...)` | the static factory — §3. |

`Options` (`prst.h:10-48`) is a bag of `std::optional`s. The fields the tests in this doc consult: `ForceFermat`, `OrderA` (the base for `-order`), `Check`/`CheckNear` (round-off error checking), `CheckStrong`/`StrongCount`/`StrongL`/`StrongL2` (Gerbicz / Gerbicz-Li), `FermatBase`, `AllFactors`, and the `Proof*`-related fields covered in `proof-system.md`.

## 3. Annotated `Run::create` — the dispatcher

`prst.cpp:463-529`, verbatim with inline notes. This is the single source of truth for "which class handles which form." It reads in code order; the first matching branch wins.

```cpp
Run* Run::create(InputNum& input, Options& options, Logging& logging, Proof* proof)
{
    // -order: only for fully-factored K*B^N+1 primes
    if (options.OrderA && !options.OrderA->empty() && options.OrderA->value() > 1)
    {
        if (input.type() != InputNum::KBNC || input.c() != 1 || !input.cofactor().empty())
        {
            logging.error("Order can be computed only for fully factored K*B^N+1 primes.\n");
            return nullptr;
        }
        return new Order(input, options, logging);
    }

    // -divides: same gate as -order (see §4.7)
    if (!options.Divides.empty())
    {
        if (input.type() != InputNum::KBNC || input.c() != 1 || !input.cofactor().empty())
        {
            logging.error("Fermat divisibility can be tested only for fully factored K*B^N+1 primes.\n");
            return nullptr;
        }
        return new FermatDivisor(input, options, logging);
    }

    // -proof cert: reuse the pre-built Proof instance (see proof-system.md)
    if (proof != nullptr && proof->op() == Proof::CERT)
        return proof;

    // -fermat, or a number no deterministic test applies to: probabilistic Fermat
    if (options.ForceFermat || input.type() == InputNum::GENERIC || std::abs(input.c()) != 1)
        return new Fermat(Fermat::FERMAT, input, options, logging, proof);

    // Proth: k*2^n+1 with k < 2^n
    if (input.type() == InputNum::KBNC && input.c() == 1 && input.b() == 2 && log2(input.gk()) < input.n())
        return new Fermat(Fermat::PROTH, input, options, logging, proof);

    input.expand_factors();
    if (!input.is_half_factored())
    {
        logging.warning("Not enough factors for an available deterministic test.\n");
        return new Fermat(Fermat::AUTO, input, options, logging, proof);
    }

    // Pocklington: c == +1, half-factored
    if (input.c() == 1)
    {
        // Simpler version of the test for small number of factors.
        if (input.type() == InputNum::KBNC && input.n() > 10 && input.factors().size() < 10)
            return new Pocklington(input, options, logging, proof);
        return new PocklingtonGeneric(input, options, logging);

    }

    // Morrison (based on Lucas chains, 2 times slower): c == -1, half-factored
    if (input.c() == -1)
    {
        // Simpler version of the test for small number of factors.
        if (input.type() == InputNum::KBNC && input.n() > 10 && input.factors().size() < 10)
            return new Morrison(input, options, logging);
        return new MorrisonGeneric(input, options, logging);

    }

    return nullptr;
}
```

The branch order matters: `-order`, `-divides`, and `-proof cert` short-circuit before any form analysis; `ForceFermat`/`GENERIC`/`|c|≠1` route to plain Fermat before the Proth check; only after `expand_factors()` + `is_half_factored()` does the dispatcher commit to a deterministic Pocklington/Morrison test. **The small-factor predicate is identical on both deterministic branches** (`prst.cpp:512` and `:522`):

```cpp
input.type() == InputNum::KBNC && input.n() > 10 && input.factors().size() < 10
```

True → the lean `Pocklington` / `Morrison`; false → the `FactorTree`-driven `PocklingtonGeneric` / `MorrisonGeneric`. The final `return nullptr` is effectively unreachable: the `|c| ≠ 1` guard at `:494` already routed everything except `c == ±1` to Fermat, so on the deterministic path either the `c==1` (Pocklington) or `c==-1` (Morrison) branch always fires. It's a defensive backstop — a `nullptr` from `create` is a hard failure in `main()` (the `-order`/`-divides` gates return it deliberately, after logging an error).

## 4. The leaf tests, one by one

Each subsection: what the constructor builds, what `run()` does, and the result line(s) it can emit.

### 4.1 `Fermat` — FERMAT / PROTH / AUTO (and the POCKLINGTON base)

**Header** (`fermat.h:11-43`): four `static const int` mode tags — `AUTO=0`, `FERMAT=1`, `PROTH=2`, `POCKLINGTON=3` — plus `int _type`, `int _a` (the base), and four task slots: `_task` (the main `MultipointExp`) and three optional `CarefulExp` helpers `_task_tail_simple` / `_task_ak_simple` / `_task_fermat_simple`.

**Constructor** (`fermat.cpp:106-308`). The shape of the exponent depends on the mode:
- A `smooth` flag (`fermat.cpp:108`) is set for `k*2^n+1` with small `k` — it lets the test exponentiate by `b=2` raised to `n` directly rather than by the full `N-1`.
- For non-FERMAT modes it walks `input.factors()` to build `exp` (and, for POCKLINGTON, a separate `exp_pocklington` / `exp_fermat`), splitting off the power of 2 into `n` (`fermat.cpp:124-145`). If there is no factor of 2 the number is even → `_a = -2` and an immediate "divisible by 2" result (`fermat.cpp:146-152`).
- It then decides the *effective* `_type`: if `log2(exp) < n` it's actually a Proth form (`fermat.cpp:160-165`); else POCKLINGTON or FERMAT.
- The base `_a` is either `genProthBase(...)` (for Proth/Pocklington, which can also detect a small algebraic factor and return negative → "divisible by %d" result, `fermat.cpp:181-187`) or `options.FermatBase` defaulting to `3` (`fermat.cpp:190, 205`).
- Finally it picks the main task class by mode and check level (`fermat.cpp:215-299`): `FastExp` (no strong check), `GerbiczCheckExp`/`FastLiCheckExp` (strong check), or `MultipointExp`/`StrongCheckMultipointExp` (when a `proof` is present, with `proof->on_point` bound as the point callback). Each registers its cost via `logging.progress().add_stage(...)`.

**`run()`** (`fermat.cpp:310-446`). The `Run::run` override (`:310`) just calls the `Proof*` overload (`:315`) with `nullptr`. The overload:
1. `_a < 0` → return immediately (a divisor was already found and reported in the ctor) (`fermat.cpp:321-322`).
2. Logs the test banner (`Proth test of …` / `Fermat probabilistic test of …`) (`fermat.cpp:324-327`).
3. Computes the `tail` correction via `_task_tail_simple` if present; a `NoInverseException` here means `N` shares a factor with the base → "is not prime, divisible by %s" (`fermat.cpp:335-363`).
4. Initializes the main task (`init_smooth` / `init_small`, or the `StrongCheck` variants), wires the proof state if present, and seeds the smooth `a^k` start value (`fermat.cpp:365-396`).
5. `_task->run()` (`fermat.cpp:400`) — the heavy exponentiation.
6. Determines `_success`: when `_task_fermat_simple` is present — the Proth/Pocklington path, where `exp_fermat != 1` (`fermat.cpp:192-193`) — its reduced result must equal `1`; otherwise `_task->result()` itself must equal `1` (`fermat.cpp:402-410`).
7. Outcome (`fermat.cpp:412-436`). Note PROTH temporarily does `*result += 1` before these comparisons and `-= 1` after (`:413, :439`):
   - PROTH with result `0` or `N` → `_prime = true`, **`"%s is prime! Time: %.1f s.\n"`**.
   - PROTH (other) or `!_success` → **`"%s is not prime. RES64: %s, time: %.1f s.\n"`**.
   - `!_prime && _success` → **`"%s is a probable prime. Time: %.1f s.\n"`**. `logging.result`'s first parameter is `bool success` (`logging.h:70`); the value passed here is `type() != PROTH && type() != POCKLINGTON` — i.e. reported as a success (**true**) for a *plain Fermat* probable prime, but **false** for a Proth/Pocklington base, since those forms either resolve primality in the branches above or hand off to the `Pocklington` subclass for the factor checks.
8. If a proof is attached, `proof->run(...)` runs the post-Fermat SAVE/BUILD (`fermat.cpp:441-442`; see `proof-system.md`).

### 4.2 `Pocklington` — small-factor, inherits `Fermat`

**Header** (`pocklington.h:6-26`): `Pocklington : public Fermat`, overriding the `Proof*` `run` overload. Adds a nested `FactorTask {int index; unique_ptr<CarefulExp> taskFactor, taskCheck;}`, a `vector<FactorTask> _tasks`, and `Giant _done` (running product of confirmed factors).

**Constructor** (`pocklington.cpp:17-45`): chains to `Fermat(Fermat::POCKLINGTON, ...)`. If the parent decided this isn't really POCKLINGTON or found a divisor (`type() != POCKLINGTON || _a < 0`) it returns early and behaves as a plain Fermat. Otherwise it builds one `FactorTask` per not-yet-done factor: `taskFactor` raises to `exp_fermat / factor`, `taskCheck` raises that by `factor` (to confirm the result is a root) (`pocklington.cpp:31-44`). Factors already marked done in the progress file are folded into `_done` instead (resume support).

**`run()`** (`pocklington.cpp:47-181`):
1. If not POCKLINGTON-mode / divisor found → delegate straight to `Fermat::run` and return (`pocklington.cpp:49-53`). This is the "Pocklington reuses Fermat::run" hand-off — it's also called for the *main stage* below even in the normal path.
2. Otherwise: open per-base child checkpoint files, log the banner, and run `Fermat::run(...)` as the probable-prime main stage (`pocklington.cpp:55-60`).
3. Loop over `_tasks` (`pocklington.cpp:63-171`): bail if `!success()`. For each factor, run `taskFactor` then `taskCheck`; a `taskCheck` result `≠ 1` means an arithmetic error → `continue` (retry). When `a^(exp/p) - 1` is nonzero, add it to the GCD batch `G`, fold the factor into `_done`, and erase the task (`pocklington.cpp:72-103`).
4. If `G` is nonempty, multiply it (via `Product`), GCD against `N`; a nontrivial GCD → **`"%s is not prime. Factor RES64: %s.\n"`** (`pocklington.cpp:105-129`). Otherwise mark those factors `done` in the progress file.
5. Halting test (`pocklington.cpp:134-139`): if all tasks consumed, or `_done²  > N` (the BLS "more than half factored" condition), set `_prime = true`.
6. If more factors are still needed, advance to the next prime base `_a`, re-setup GWnum, and re-run `Fermat::run` with the new base (`pocklington.cpp:141-170`). **If a proof is active this restart is impossible** → "Pocklington test needs to restart, disable proofs to proceed." + abort (`pocklington.cpp:142-146`).
7. On `_prime`, emit **`"%s is prime!\n"`** — note: *no* `Time:` suffix here, unlike every other prime line (`pocklington.cpp:173-177`).

### 4.3 `PocklingtonGeneric` — FactorTree-driven

**Header** (`pocklington.h:30-47`): `Run` subclass. Holds a `SubLogging _logging`, `Giant _done`, `std::set<int> _done_factors`, a `unique_ptr<FactorTree> _tree`, and `int _a`. Used when the candidate has many factors or isn't plain `KBNC` (e.g. `n!+1`, `n#+1`).

**Constructor** (`pocklington.cpp:183-254`): recovers `_done_factors` from the comma-separated `"factors"` progress param (resume), accumulates the product of done factors into `_done` and the remaining exponent into `exp` (chunking giant multiplies at 8192 limbs), then builds the `SubLogging` with **`_logging->progress().set_parent(nullptr)`** and calls `create_tasks`. The `set_parent(nullptr)` detaches the inner progress from the parent — the parent's stage is advanced manually via `logging.progress().update(...)` in `run()`. This is the time-accounting seam; see §7.

> **`create_tasks` is almost entirely commented out** (`pocklington.cpp:256-325`). Only the live prologue (`:260-271`) matters: it groups the not-yet-done factors into leaf `FactorTree` nodes and builds the balanced `_tree`, then assigns `_tree->exp()`. The real task construction does **not** happen here — it is inlined into the `run()` stack walk (`:445-461`). Don't be misled by the large commented block describing a per-node `task()`; it's dead reference code.

**`run()`** (`pocklington.cpp:327-612`): a depth-first walk of `_tree` over an explicit `stack`/`stack_value`, restarting with a new base on each outer-loop pass until `_prime` or a factor is found (`:358`):
- For each node with a nontrivial `exp`, pick the task by bitlen and stack depth: `CarefulExp` (<32 bits), else `FastLiCheckExp`/`FastExp` at the root, else `LiCheckExp`/`SlidingWindowExp` deeper (`:446-461`). Strong-check vs. plain is chosen by `_options.CheckStrong` (default true). Run it under `_logging`, push the result onto `stack_value` (`:495-498`).
- At a leaf factor: if the residue `≠ 1` and not yet `_success` → **`"%s is not prime. RES64: %s, time: %.1f s.\n"`** (`:512`); if `≠ 1` after success it's an arithmetic error → restart; else mark `_success`, and if the popped parent value `≠ 1` add it to the GCD batch `G` (`:506-559`). Batches of 100 are flushed via the `test_G` lambda.
- `test_G` (`:374-426`) products the batch, GCDs against `N` → **`"%s is not prime. Factor RES64: %s.\n"`** on a hit (`:397`), else records the factors done and checks the `_done² > N` halting condition → `_prime`.
- On a clean pass with factors still outstanding, advance `_a` to the next prime (skipping bases with the wrong Kronecker symbol for factor 0), re-setup GWnum, rebuild tasks, and loop (`:579-598`).
- On `_prime`: **`"%s is prime! Time: %.1f s.\n"`** (`:606`).

### 4.4 `Morrison` — small-factor (`c == -1`)

**Header** (`morrison.h:12-36`): `Run` subclass for `k*b^n-1`, `n!-1`, `n#-1`. The math is Lucas-sequence based (BLS Theorem 14; see the comment block at `morrison.cpp:14-18`) and ~2× slower than Pocklington. Holds `_task` (a `LucasVMul` — concretely `LucasVMulFast` without strong check or `LucasUVMulFast` with), `_taskCheck` (`LucasVMulFast`), a `vector<FactorTask>` (each with `LucasVMulFast taskFactor/taskCheck`), and the Lucas params `int _P`, `bool _negQ`.

**Constructor** (`morrison.cpp:20-169`): walks factors to build the Lucas exponent. The power of 2 determines `_negQ` (`Q=-1` iff `2^n` with `n>1`; for `Q=-1`, factor 2 is "tested for free") (`morrison.cpp:43-49`). If `log2(exp) < n` it's the LLR special case → `_name = "Morrison (LLR) test"`, factor tasks cleared (`:86-92`). With `CheckStrong` (default) it builds a `LucasUVMulFast` with Gerbicz checks; otherwise a plain `LucasVMulFast` (`:97-109`). It then populates the Lucas chains with the per-prime multipliers (`:119-165`) and picks the smallest `_P` whose discriminant has Kronecker symbol ≠ 1 (`:167`). Even number (no factor of 2) → "divisible by 2" + early return (`:70-76`).

**`run()`** (`morrison.cpp:171-355`): outer loop over candidate `_P` values until `_prime` (`:186`):
1. On a restart, advance `_P` to the next valid one and re-stage (`:188-202`).
2. GCD `4·P·(P²∓4)` against `N`; nontrivial → **`"%s is not prime. Factor RES64: %s, time: %.1f s.\n"`** (`:210-213`).
3. Log banner, init and run the main Lucas task; if `_taskCheck` is present, chain it to reduce the V-value (`:217-242`).
4. The reduced V must equal `2` (for `Q=1`) or `0` (for `Q=-1`); otherwise → **`"%s is not a probable prime. Have you run Fermat test first? RES64: %s, time: %.1f s.\n"`** (`:245-251`). Set `_success`.
5. If there are factor tasks, compute per-factor V-values (single factor: divide the main V; multiple: run each `taskFactor`/`taskCheck` chain), accumulate the differences into a GCD batch, and check the half-factored / `square(done) > N` condition; a nontrivial GCD → **`"%s is not prime. Factor RES64: %s, time: %.1f s.\n"`** (`:255-340`).
6. `_prime = true` → **`"%s is prime! Time: %.1f s.\n"`** (`:341-348`).

### 4.5 `MorrisonGeneric` — FactorTree-driven (`c == -1`)

**Header** (`morrison.h:38-53`): mirrors `PocklingtonGeneric` (SubLogging, `_done`, `_done_factors`, `_tree`) plus `std::vector<int> _dac_index` (precomputed Lucas-chain "DAC" addition-chain indices per small prime) and `_P`/`_negQ`.

**Constructor** (`morrison.cpp:357-458`): **forces `CheckStrong = true`** — without the strong check this form is unsupported and it logs an error (`:359-362`). Recovers `_done_factors`, precomputes `_dac_index` for each small prime factor (`:403-419`), builds the `FactorTree` of odd factors (factor 2 sets `_negQ`), assigns `_tree->exp()`, and sets up the detached `SubLogging` (`set_parent(nullptr)`, same as PocklingtonGeneric).

**`run()`** (`morrison.cpp:460-857`): structurally the same stack walk as `PocklingtonGeneric::run`, but over `LucasVMulFast::State` stack values and with two extra wrinkles:
- **On-disk stack resume.** Before the walk it reads any persisted `stack_value` States from child checkpoint files, and the `write_values` lambda (`:579-619`) serializes the whole stack to disk on the `DISK_WRITE_TIME` cadence and on abort. (Pocklington's generic path doesn't persist its stack this way.)
- Task choice per node: `LucasVMulFast` (<32 bits, using `_dac_index` when available), `LucasUVMulFast` at the root, `LucasUVMul` deeper (`:662-681`).
- Leaf condition: V must equal `0` (`Q=-1`) or `2` (`Q=1`); mismatch before success → **`"%s is not a probable prime. Have you run Fermat test first? RES64: %s, time: %.1f s.\n"`** (`:751`); the GCD-batch hit in `test_G` → **`"%s is not prime. Factor RES64: %s.\n"`** (`:553`).
- On `_prime`: **`"%s is prime! Time: %.1f s.\n"`** (`:851`).

### 4.6 `Order` — multiplicative order

**Header** (`order.h:10-47`): `Run` subclass invoked by `-order` for fully-factored `K·B^N+1` primes. Computes `ord(a) mod N`. Holds `_factors` (remaining factors of `N-1`), `_order` (computed prime-power orders), a main `_task` (`MultipointExp`), `_tasks_smooth` (per-prime smooth exponentiations), `_task_check` (a `CarefulExp` verifying `a^(N-1) == 1`), per-factor `_tasks` (`FactorTask{Giant b; int ord, n; CarefulExp task_sub, task_factor;}`), and `int _sub = 30` (the small-power cutoff). A protected named constructor and the extracted `create_smooth_task` (`order.cpp:25-50`) plus a `result()` accessor exist for the `FermatDivisor` subclass (§4.7).

**Constructor / `create_tasks`** (`order.cpp:13-119`): `create_tasks` splits the factorization into a "smooth" part (prime powers with exponent `> _sub`, or factor 2, get their own `GerbiczCheckExp`/`MultipointExp` via `create_smooth_task`, with an `on_point` callback) and a combined "non-smooth" part that goes into the single main `_task` (`FastLiCheckExp`/`LiCheckExp`/`FastExp`/`SlidingWindowExp` by check level and base size) (`order.cpp:52-103`). It then builds the per-factor `task_sub`/`task_factor` `CarefulExp`s that peel one prime at a time (`:105-119`). `on_point` (`:121-129`) throws `TaskAbortException` the moment a smooth exponentiation hits `1`, recording the breakpoint index.

**`run()`** (`order.cpp:131-302`): loop while `_factors` non-empty (`:139`):
1. Run the main `_task` to raise `a` to the non-smooth exponent (`:144-159`).
2. Run each smooth task; reaching `1` removes that factor, else the `on_point` abort captures the exact power where it became `1` (`:165-199`).
3. If everything collapsed to `1`, halve the remaining powers by `_sub` and rebuild (`:201-208`).
4. Verify `a^(N-1) == 1` via `_task_check`; failure → `"%s is not prime."` + abort (the candidate wasn't actually prime) (`:212-220`).
5. For each factor, peel primes via `task_sub`/`task_factor` until the value stops being `1`, recording the order `ord` (`:222-258`).
6. After the loop, assemble the order as a factor product, optionally collapsing to `(N-1)/divisor`, and emit the result: **`"ord(%s) mod %s = %s.\n"`** (`order.cpp:297`). Checkpoints use the `.ord<fingerprint>` filename suffix (`prst.cpp:341-342`).

### 4.7 `FermatDivisor` — Fermat/GF/xGF divisibility search (`-divides`)

**Header** (`order.h:49-80`): `FermatDivisor : public Order` — reuses Order's smooth-task machinery. A nested `Base {base, str, power, task, sub_val, exp, val}` holds per-base state.

**CLI**: `-divides {f | gf | xgf} [limit]` (`prst.cpp:147-163`; also available under `-batch`, `batch.cpp:78-94`). `f` searches Fermat numbers only (base 2); `gf`/`xgf` search bases 2..`limit` (default 12). The dispatcher gate is the same as `-order`'s (`prst.cpp:478-486`); the constructor additionally asserts the input's first factor is 2 (`order.cpp:306`). Checkpoints use the `.div` filename suffix (`prst.cpp:343-344`) with per-base/per-power child files (`order.cpp:346-347`); resume goes through `base_<b>` progress params. `base.task->write_state()` is called right after `run()` (`order.cpp:361`) — a live caller of the public `Task::write_state()`.

**`run()`** (`order.cpp:386-542`): for each prime base `b`, raise `b` to `2^(n−30)` with a checkpointed smooth task (`create_smooth_task`) and finish the last 30 squarings with a `CarefulExp`; composite bases are assembled from their prime factors' values with `Product::mul` instead of running their own exponentiation. If `b^(2^m) ≡ −1 (mod P)` for some `m` (found by squaring toward `N−1`), P divides `F(m)` (base 2) or `GF(m,b)` — reported only when `perfect_power(b) == 1`, since perfect-power bases are redundant. In `xgf` mode, pairs `(a,b)` with `gcd(a,b) = 1` (and not both `perfect_power` even) are checked for `a^(2^m) + b^(2^m) ≡ 0 (mod P)`; misaligned powers trigger a fix-up exponentiation on the base `a/b (mod P)` (`order.cpp:495-502`).

**Result lines** (`order.cpp:528-538`): **`"%s divides %s.\n"`** with a list like `F(2543548), GF(2543549,3), xGF(2543549,3,2)`, or **`"%s no divisible numbers found.\n"`**. Not a primality verdict — like `Order`, this mode assumes the input is already known prime.

## 5. `FactorTree` — construction and walk

`FactorTree` (`pocklington.h:49-108`) is the structure both `*Generic` tests walk. It's a binary tree whose leaves are individual factors (each carrying `_index` into `input.factors()`) and whose internal nodes carry the *product* of their subtree's exponents. Walking it lets the test raise the base to progressively larger partial exponents, reusing each parent's result as the starting point for its children — far cheaper than recomputing `a^(exp/p)` per factor independently.

The interesting constructor is the vector overload (`pocklington.h:58-94`), which does a **balanced bottom-up merge**:

```cpp
FactorTree(std::vector<std::unique_ptr<FactorTree>>& tree) : _index(-1)
{
    if (tree.size() == 1)
    {
        _left = std::move(tree[0]);
        return;
    }
    for (auto it = tree.begin(); it != tree.end(); it++)
        (*it)->_left.reset(new FactorTree((*it)->exp(), (*it)->index()));

    while (true)
    {
        std::vector<std::unique_ptr<FactorTree>> next;
        for (auto it = tree.begin(); it != tree.end(); )
        {
            auto& a = *it;  it++;
            if (it == tree.end())
                next.push_back(std::move(a));          // odd one out, carried up
            else
            {
                auto& b = *it;  it++;
                swap(a->exp(), b->exp());               // children hold each other's exp
                std::swap(a->_index, b->_index);
                if (tree.size() == 2)
                {
                    _left = std::move(tree[0]);
                    _right = std::move(tree[1]);
                    return;
                }
                next.emplace_back(new FactorTree(a->exp()*b->exp()));  // internal node = product
                next.back()->_left = std::move(a);
                next.back()->_right = std::move(b);
            }
        }
        tree = std::move(next);
    }
}
```

Accessors (`pocklington.h:96-101`): `index()` (leaf factor index, `-1` for internal nodes), `exp()` (the node's exponent), `is_factor()` (`!_left && !_right`, a leaf), `is_last()` (`!_right`, no right child), and `left()`/`right()`.

The walk (in both `*Generic` `run()`s) is an explicit-stack DFS:

```cpp
stack.push_back(_tree.get());
while (!stack.empty())
{
    // ... if node has a nontrivial exp: build+run a task, push result onto stack_value ...
    if (!stack.back()->is_factor() && !stack.back()->left()->exp().empty())
        stack.push_back(stack.back()->left());
    else if (!stack.back()->is_last() && !stack.back()->right()->exp().empty())
        stack.push_back(stack.back()->right());
    else
    {
        stack.pop_back();
        if (!stack_value.empty())
            stack_value.pop_back();
    }
}
```

`stack_value` shadows `stack`: each entry is the partial exponentiation result feeding that node's children. PocklingtonGeneric pushes `Giant`s; MorrisonGeneric pushes `LucasVMulFast::State`s and additionally persists/restores them across restarts.

## 6. Common patterns

- **Small-factor vs. `*Generic` is one predicate.** `KBNC && n>10 && factors.size()<10` (`prst.cpp:512, 522`). The lean classes work directly on the parent `Logging`; the `*Generic` classes detach a `SubLogging` and walk a `FactorTree`. When debugging "which code runs for this input," evaluate that predicate first.
- **Pocklington reuses Fermat by inheritance.** `Pocklington : Fermat` calls `Fermat::run` for the probable-prime stage, then layers per-factor Pocklington checks on top. None of the others inherit a sibling test.
- **GCD-batch + restart-with-new-base.** Pocklington (small and generic) and both Morrison variants share the shape: accumulate `result-1` (or the Lucas difference) into a vector `G`, `Product` them, GCD against `N`, and if the candidate isn't yet "more than half proved," advance the base/`_P` to the next prime and re-run. The halting condition is always a variant of `square(_done) > N`.
- **Resume via the progress file.** The `*Generic` classes serialize completed factor indices into the `"factors"` progress param and skip them on construction; small Pocklington uses per-factor `"factor<i>" = "done"` params. `_a`/`_P` are likewise stashed as params for resume.
- **`add_stage` / `next_stage` / `update`.** Every task's cost is registered with `logging.progress().add_stage(task->cost())` before it runs, and `progress().next_stage()` advances after. The `*Generic` tests reset `_logging->progress() = Progress()` per node and push progress to the parent manually with `logging.progress().update(...)`.

## 7. Pitfalls

- **Result-line wording is a user-visible contract.** mersenneforum users and BOINC log parsers match these exact strings (`is prime!`, `is not prime`, `RES64:`, `time: %.1f s.`). Note the *inconsistencies already in the tree*: small `Pocklington` prints `"%s is prime!\n"` with **no** `Time:` suffix (`pocklington.cpp:175`), while every other prime line includes `Time: %.1f s.`. Don't "fix" wording casually — see lay-of-the-land note 3.
- **`*Generic` time accounting is subtle.** Both generic tests detach their `SubLogging` progress (`_logging->progress().set_parent(nullptr)`, `pocklington.cpp:250`, `morrison.cpp:456`) and feed the parent's stage manually. Generic time accounting is subtle; the fix for an earlier mis-reported-time bug landed in `d14b2f7`. Read `update`/`next_stage` carefully before changing this.
- **`create_tasks` in PocklingtonGeneric is mostly dead code.** The live work is the small prologue (`pocklington.cpp:260-271`); the actual task objects are built inline in `run()`'s stack walk. Editing the commented block does nothing.
- **Proofs only work on the Fermat family.** `Pocklington::run` aborts with "disable proofs to proceed" the moment it needs a base restart (`pocklington.cpp:142-146`), and the `*Generic`/`Morrison`/`Order` constructors take no `Proof*`. A new deterministic test that wants proof support has to route through the Fermat/`MultipointExp` machinery.
- **`MorrisonGeneric` silently forces the strong check.** It overrides `CheckStrong` to `true` regardless of the option (`morrison.cpp:362, 472`). A user passing `-CheckStrong 0` for a generic Morrison candidate gets strong-checked anyway (with an error logged at construction).
- **The small-factor predicate can route a "Proth-looking" input to Pocklington.** Proth is checked *before* `expand_factors()`; only `k*2^n+1` with `b==2` and `log2(k)<n` is Proth. A `k*b^n+1` with `b≠2`, or `k≥2^n`, falls through to the Pocklington branch even though a user might expect "Proth." Trace `Run::create` top-to-bottom.

## 8. Quick reference

| Class | Base | Dispatch condition (`Run::create`) | Main task type(s) | Prime line | Not-prime / factor line |
|---|---|---|---|---|---|
| `Fermat` (FERMAT) | `Run` | `ForceFermat` / `GENERIC` / `\|c\|≠1` | `FastExp` / `FastLiCheckExp` / `MultipointExp` | `is a probable prime. Time:` | `is not prime. RES64: …, time:` |
| `Fermat` (PROTH) | `Run` | `KBNC, c=1, b=2, log2(k)<n` | same | `is prime! Time:` | `is not prime. RES64: …, time:` |
| `Pocklington` | `Fermat` | `c=1`, half-factored, `KBNC, n>10, factors<10` | `Fermat::run` + per-factor `CarefulExp` | `is prime!` *(no time)* | `is not prime. Factor RES64:` |
| `PocklingtonGeneric` | `Run` | `c=1`, half-factored, else | tree of `CarefulExp`/`Li`/`SlidingWindow` | `is prime! Time:` | `is not prime. RES64: …, time:` / `Factor RES64:` |
| `Morrison` | `Run` | `c=-1`, half-factored, `KBNC, n>10, factors<10` | `LucasVMulFast` / `LucasUVMulFast` | `is prime! Time:` | `is not prime. Factor RES64: …, time:` / `is not a probable prime. …` |
| `MorrisonGeneric` | `Run` | `c=-1`, half-factored, else | tree of `LucasVMulFast`/`LucasUVMul(Fast)` | `is prime! Time:` | `is not prime. Factor RES64:` / `is not a probable prime. …` |
| `Order` | `Run` | `-order` on fully-factored `KBNC, c=1` | `MultipointExp` + `CarefulExp` | `ord(a) mod N = …` | `is not prime.` (error) |
| `FermatDivisor` | `Order` | `-divides` on fully-factored `K*2^N+1` | smooth `MultipointExp`/`GerbiczCheckExp` + `CarefulExp` + `Product` | `divides F(…), GF(…), xGF(…).` | `no divisible numbers found.` |

## 9. Open questions / non-coverage

- **The `_smooth` exponentiation path.** `Fermat`'s ctor toggles a `smooth` mode for small-`k` `k*2^n+1`, changing the exponent representation and which `init_*` is called. The mechanics are in `exponentiation-algorithms.md`; this doc only notes where the flag is read.
- **The math.** BLS Theorem 14 (Morrison), the Pocklington half-factored condition, Gerbicz-Li strong checks, and the Lucas-chain / DAC addition-chain construction are described operationally here but not derived. See `framework/docs/mult_*.pdf` and eprint 2023/195; deferred indefinitely (Tier 3 in lay-of-the-land) unless modifying the algorithms.
- **Multi-`Task` parallelism.** Every `run()` here is single-threaded above GWnum's FFT layer. Parallelizing the `*Generic` factor walk would require locking the entire `Logging`/`Progress` subsystem — see `logging-and-progress.md` §12 and `task-lifecycle.md` Pitfall E.
- **Task internals.** `init_*`, `commit_*`, the restart/checkpoint loop, and `on_state` cadence are `task-lifecycle.md`'s domain; this doc treats `task->run()` as a black box.
