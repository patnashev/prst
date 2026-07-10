# Test harness — deep dive

`PRST -test <subset>` is the built-in regression suite. It re-runs known primality tests against **precomputed expected answers** — a 64-bit residue (`res64`) and a 64-bit certificate hash (`cert64`) per candidate — and aborts the moment one disagrees. So it's two things at once: a correctness check on the test algorithms *and* the project's own machine-readable definition of "the right answer." If you change anything in the Fermat/Proth/Pocklington/Morrison or proof code, this is what tells you whether you broke it.

It also bundles two non-primality suites: **`RootsTest`** (the `error` subset) — which *forges* proof points to confirm the proof system rejects tampering — and **`ABCParserTest`** (the `abc_parser` subset) — unit tests for the batch-file parser. Separately, `framework/arithmetic/test.cpp` is a **standalone `main()`** smoke test of the arithmetic layer; it is *not* wired into `-test` (it's its own executable).

`testing_main` is an alternate entry point like `-batch`/`-boinc`/`-net`: `prst.cpp:197` dispatches `-test` into it via `exit(testing_main(...))`.

Source files:
- `src/testing.h`, `src/testing.cpp` (`testing_main`, the `Test`/`DeterministicTest` classes, `RootsTest`, `ABCParserTest`)
- `src/test.data` (a 5,560-line `#include`d C source file — the canned vector arrays)
- `framework/arithmetic/test.cpp` (the standalone arithmetic smoke test)

Prereqs / companions: `proof-system.md` (`Test::run` exercises the full SAVE→BUILD→CERT pipeline; `RootsTest` attacks the roots-of-unity defense from its §7), `run-hierarchy.md` (`DeterministicTest` goes through `Run::create`), `inputnum-parsing.md` (every `Test` wraps an `InputNum`), `abc-batch.md` (`ABCParserTest` drives `parse_batch_file`), `config-dsl.md` (the option chain), `arithmetic-foundation.md` (`arithmetic/test.cpp` exercises `Giant`/`GWNum`/`ReliableGWArithmetic`/Edwards).

## 1. Data model

A **`Test`** (`testing.h:37`) wraps an `InputNum` plus its expected `res64`/`cert64`. It's built from one of four POD vector structs, each a row in `test.data`:

| Struct | Fields | Candidate form |
|---|---|---|
| `NTest` | `n, res64, cert64` | a fixed `k`/`b`/`c` with varying `n` (e.g. `3*2^n+1`) |
| `BTest` | `b, res64, cert64` | varying base `b` at fixed `n` (generalized Fermat `b^8192+1`) |
| `KBNCTest` | `sk, sb, n, c, res64, cert64` | an explicit `k*b^n+c` (strings for `k`,`b`) |
| `FreeFormTest` | `s, bitlen, res64, cert64` | an arbitrary expression string |

`DeterministicTest : Test` (`testing.h:93`) is the variant for the deterministic tests (Pocklington/Morrison); it overrides `display_text()` to label which test will run and `run()` to drive `Run::create` instead of a Fermat+proof. `TestLogging` (`testing.h:119`) is a `Logging` whose `heartbeat()` emits a progress line every `PROGRESS_TIME` seconds.

The vector arrays in `test.data` (each terminated by a zero sentinel):

| Array | Built as | Family |
|---|---|---|
| `Test321Plus` / `Test321Minus` | `Test(3, 2, n, ±1)` | `3*2^n±1` ("321") |
| `TestBase5Plus` / `TestBase5Minus` | `Test(2, 5, n, ±1)` | `2*5^n±1` |
| `TestGFN13` / `TestGFN13More` | `Test(1, b, 8192, 1)` | generalized Fermat `b^(2^13)+1` |
| `TestSpecial` / `TestPrime` | `Test(KBNCTest)` | assorted `k*b^n+c` |
| `TestFreeForm` | `Test(FreeFormTest)` | arbitrary expressions |
| `Test109208Base5Plus` / `Test100186Base5Minus` | `Test(k, 5, n, ±1)` | the `slow` large-`k` runs |

`res64`/`cert64 == 0` in a row means **"don't verify, just record"** — `Test::run` fills them from the run instead of comparing (used for the ad-hoc custom-expression path). A `FreeFormTest` with `res64 == 0` is special-cased in the ctor to `res64 = 1` (`testing.h:68`), meaning "expected prime."

## 2. `Test::run` — the central method

`testing.cpp:293-399`. This is a full proof round-trip with a cross-check after every phase — the densest validation in the codebase. Annotated skeleton:

```cpp
void Test::run(Logging& logging, Options& global_options, GWState& global_state)
{
    Options options;                       // fresh, with Check* inherited
    options.ProofSecuritySeed = "12345";
    options.RootOfUnityCheck = false;
    int proof_count = 16;

    Proof proof(Proof::SAVE, proof_count, input, options, file_cert, logging);
    Fermat fermat(Fermat::AUTO, input, options, logging, &proof);
    GWState gwstate; gwstate.copy(global_state); input.setup(gwstate);

    // --- phase 1: SAVE (run the Fermat test, recording proof points) ---
    proof.init_files(&file_proofpoint, &file_proofproduct, &file_cert);
    fermat.run(gwstate, file_checkpoint, file_recoverypoint, logging, &proof);
    if (res64 == 0)  res64  = fermat.success() ? 1 : std::stoull(fermat.res64(), …);   // record mode
    if (cert64 == 0) cert64 = std::stoull(proof.res64(), …);
    if (fermat.success() != (res64 == 1))                 throw …;  // primality mismatch
    if (!fermat.success() && stoull(fermat.res64()) != res64) throw …;  // RES64 mismatch
    if (stoull(proof.res64()) != cert64)                  throw …;  // raw certificate mismatch

    // --- phase 2: BUILD (compress points into a certificate) AT A DIFFERENT FFT ---
    gwstate.done(); gwstate.next_fft_count = 1; input.setup(gwstate);
    Proof proof_build(Proof::BUILD, proof_count, input, options, file_cert, logging);
    proof_build.calc_points(…); proof_build.init_files(…);
    fermat.run(gwstate, file_checkpoint, file_recoverypoint, logging, &proof_build);
    if (… build raw cert … != cert64)  throw …;            // build raw certificate mismatch
    if (fermat.success() != (res64 == 1)) throw …;          // + primality / RES64 re-checks

    // --- phase 3: CERT (re-verify the certificate) BACK AT THE DEFAULT FFT ---
    gwstate.done(); gwstate.next_fft_count = 0; input.setup(gwstate);
    Proof proof_cert(Proof::CERT, 0, input, options, file_cert, logging);
    proof_cert.run(gwstate, file_checkpoint, file_recoverypoint, logging);
    if (proof_cert.res64() != proof_build.res64())  throw …; // certificate mismatch
}
```

Seven comparison-and-abort checks across the three phases — primality, RES64, and raw-certificate after SAVE (`:337,:342,:347`); the same three after the BUILD re-run (`:364,:369,:374`); and the final SAVE-vs-BUILD certificate match (`:387`) — and the BUILD phase deliberately runs **at a larger FFT** (`next_fft_count = 1`) than SAVE/CERT, so a passing test also confirms the result is **FFT-independent** (the same number computed two different ways must give the same residue and certificate). A mismatch anywhere throws `TaskAbortException`, which `testing_main` catches and reports as `Failed test: <text>`. (This re-runs the Fermat from scratch at the second FFT; it validates *result* consistency, not checkpoint portability — see `arithmetic-foundation.md` §9.)

`DeterministicTest::run` (`testing.cpp:401-458`) is the simpler analogue: `expand_factors()`, require `is_half_factored()` (else `runtime_error("Not enough factors.")`), `Run::create` → `run` → check `prime() == (res64 == 1)` and the RES64.

## 3. `testing_main` and the subsets

`testing_main` (`testing.cpp:30-291`) parses options with the usual `Config` chain (see `config-dsl.md`) and takes the subset name via `default_code`. It builds a `list` of `(name, SubLogging, deque<Test>)`, registers each test's `cost()` (`bitlen^1.5`) as a progress stage, then runs them, each under a `SubLogging` that suppresses per-test output (down to `LEVEL_ERROR`) unless the log level is `-log debug` — note `-d` (which sets `LEVEL_INFO`) does **not** unsuppress it, only `-log debug` does (`testing.cpp:269`). On success: `All tests completed successfully.`; on any mismatch: `Failed test: <display_text>` and `PRST_EXIT_FAILURE`.

| Subset | Contents |
|---|---|
| `321plus` / `321minus` | `3*2^n±1` |
| `b5plus` / `b5minus` | `2*5^n±1` |
| `gfn13` | `b^8192+1` |
| `special` | `TestSpecial` (assorted) |
| `freeform` | `TestFreeForm` (each also run as a `DeterministicTest` when `res64==1`) |
| `deterministic` | `TestPrime` entries via `Run::create` (only `c=-1`, or `c=1` with `b≠2`) |
| `prime` | `TestPrime` via Fermat+proof |
| `error` | `RootsTest` (the proof-forgery defense suite) |
| `all` | `321plus + 321minus + b5plus + b5minus + gfn13 + special + error + freeform + deterministic + prime` |
| `slow` | `gfn13more + 100186b5minus + 109208b5plus` |
| `abc_parser` | `ABCParserTest`, returns immediately |
| *(anything else)* | treated as a **custom expression** — parsed and run once in record mode (`res64=0`) |

So `PRST -test "1234567*2^9999+1"` runs that one candidate through the full proof pipeline without a pre-stored answer.

## 4. The two non-primality suites

**`RootsTest`** (`testing.cpp:460-759`, the `error` subset) verifies the proof system's *forgery defenses* — it's named for testing error/attack **detection**, not for being broken. For `3*2^353+1`, `960^128+1`, and `2*5^178-1` it builds a genuine proof, then tampers:
- writing a **random** value over a proof point → the rebuilt certificate must *not* match (catches naive tampering; "Certificate failure" if it wrongly matched);
- multiplying proof points by a **root of unity** `R` (a forgery that stays internally consistent) → via the `test_attack` lambda: the Fermat must report *not prime* ("Attack failure" otherwise), the tampered certificate *does* still match the build (the forgery is self-consistent), **but** the `Proof::ROOT` roots-of-unity check must flag it (`taskRoot()->result() == 1` = attack detected; "Roots of unity check failed to detect the attack" otherwise).

This is the executable proof that the §7 defenses in `proof-system.md` actually work.

**`ABCParserTest`** (`testing.cpp:782+`, the `abc_parser` subset) is a self-contained unit test of the batch parser: it writes temp files and asserts on `detect_format` and `parse_batch_file` output — candidate counts, expanded expressions, `k_value`s, comment/blank-line skipping, ABCD delta accumulation (`ABCD $a*2^1000+1 [3]` + deltas `2 4 6` → 4 candidates starting `3*2^1000+1`), out-of-bounds returns. A `check(condition, name)` lambda tallies failures; it returns the failure count (0 = pass) rather than throwing.

**`framework/arithmetic/test.cpp`** is a different animal: a standalone program with its own `main()` (not reachable via `-test`). It's a developer smoke test that exercises the arithmetic layer directly — Edwards/Montgomery curve arithmetic and `gen_curve`, `Giant` operators and big-number string round-trips, `LucasV` sequences, a `ReliableGWArithmetic` Proth squaring with the restart loop (`224027*2^99763+1`), and `Poly`/`FFT` polynomial multiplication — printing results to `cout` for manual inspection. It's compiled separately and run by hand, not part of the regression suite.

## 5. Lifecycle

1. **Parse.** `testing_main` runs the `Config` chain; `default_code` captures the subset name.
2. **Materialize.** For the requested subset(s), walk the matching `test.data` array(s) and `emplace_back` a `Test`/`DeterministicTest` per row into a per-subset deque. (`abc_parser` short-circuits to `ABCParserTest` and returns.)
3. **Stage.** Each test's `cost()` is added as a progress stage so the overall bar is meaningful.
4. **Run.** Per subset, per test: `display_text()` (also re-parses the input and labels the test type), then `test->run(...)` under a `SubLogging`. The `error` subset runs `RootsTest` instead.
5. **Verdict.** First mismatch → `TaskAbortException` → `Failed test: <text>` + exit 1. All pass → `All tests completed successfully.` + exit 0. Ctrl-C → `Test aborted.`

## 6. Pitfalls

- **`res64`/`cert64` are the oracle — changing the algorithm means regenerating them.** They're precomputed answers baked into `test.data`. A legitimate change to how a residue or certificate is computed (not a bug) will make every affected row "fail"; the fix is to regenerate the data, not to relax the check. Treat a mass failure after an intentional change as "update the oracle," a single failure as "you broke something."
- **`0` means "record," not "expected zero."** A row with `res64 == 0` (or the custom-expression path) is *not* checked — `Test::run` fills it from the run. Only the special `FreeFormTest` ctor remaps `0 → 1` ("expect prime"). Don't read a `0` in the data as "the residue should be zero."
- **The `error` subset is the defense suite, not a list of broken cases.** It's `RootsTest`, which *intentionally* forges proofs and asserts they're caught. A failure there means a proof-system *defense* regressed, which is more serious than a normal RES64 mismatch.
- **`DeterministicTest` throws on under-factored input.** It calls `is_half_factored()` and `throw`s `runtime_error("Not enough factors.")` rather than skipping — a `test.data` row that isn't actually half-factored aborts the whole run.
- **`framework/arithmetic/test.cpp` is not `-test`.** It has its own `main()`, isn't built into the `prst` binary, and reports via `cout` for eyeballing (no automated pass/fail). Don't expect `-test` to cover it, and don't expect it to be a clean regression gate.
- **BUILD runs at a different FFT on purpose.** The `next_fft_count = 1` in phase 2 isn't a leftover — it's the cross-FFT consistency check. Removing it would weaken the suite (it would no longer catch FFT-size-dependent result bugs).

## 7. Quick reference

| You want to… | Do |
|---|---|
| Run the standard regression | `PRST -test all` |
| Run the long/large-`k` cases | `PRST -test slow` |
| Test one form's family | `PRST -test 321plus` (or `b5minus`, `gfn13`, …) |
| Test the proof-forgery defenses | `PRST -test error` |
| Unit-test the batch parser | `PRST -test abc_parser` |
| Run one ad-hoc candidate (no oracle) | `PRST -test "<expression>"` |
| Add a regression vector | append a row (with computed `res64`/`cert64`) to the right array in `test.data` |
| Find what "correct" is for a candidate | the matching `res64`/`cert64` in `test.data` |

## 8. Open questions / non-coverage

- **How `test.data` is generated.** The 5,560 rows of `res64`/`cert64` are clearly machine-produced (likely from trusted prior runs / external references like the 321 and base-5 search projects), but the generator/provenance isn't in this repo. Regenerating after an intentional algorithm change is a real task with no documented tool here.
- **The exact `RootsTest` forgery constructions.** This doc explains *what each attack checks* (random tamper rejected; root-of-unity forgery caught by the ROOT check); the specific algebra of each `R` (the `FastExp` exponents, the per-point multiplications) is reproduced from the code but its number-theoretic *why* belongs with `proof-system.md` §7 and `math-and-theorems.md`.
- **`arithmetic/test.cpp` assertions are by-eye.** It prints values rather than asserting; knowing whether its output is "right" requires the expected values (mostly implicit — e.g. `prime`, matching `==` booleans). A fuller treatment belongs with `arithmetic-foundation.md` / `curves-and-polynomials.md`, which own those classes.
- **Per-suite coverage rationale.** *Why* these specific `n`/`b` values were chosen as the regression set (FFT-boundary cases, historically-buggy sizes, etc.) isn't recorded; the data is treated as an opaque trusted oracle.
