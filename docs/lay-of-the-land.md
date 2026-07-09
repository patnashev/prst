# PRST — Lay of the Land

A high-level orientation to the PRST application for contributors new to the codebase. For deeper dives
into the test-orchestration subsystems see the per-subsystem docs in this folder. PRST is built on a
shared library, the **arithmetic framework** (`patnashev/arithmetic`, consumed here as the `framework/`
submodule); that library has its own orientation and deep-dives in its repo's `docs/` folder — see
[The framework](#the-framework-shared-library) below.

## What PRST is

PRST is a C++ command-line utility for **primality testing of very large numbers**. Written by Pavel Atnashev (`patnashev/prst` on GitHub), built on George Woltman's GWnum multiplication library. It is best at systematic searches of numbers in popular forms — Proth, Thabit, generalized Fermat, factorials, primorials, arbitrary candidates — and is typically driven by a sieving stage upstream of it. (Mersenne numbers go to GIMPS instead.)

It implements:
- **Fermat** probabilistic test (any number, not definitive)
- **Proth** test of `k*2^n+1`
- **Pocklington** test of `k*b^n+1`, `n!+1`, `n#+1`
- **Morrison** test of `k*b^n-1`, `n!-1`, `n#-1`
- **LLR** as a special case of Morrison

Plus a verification-grade proof system (`-proof save/build/cert`) that lets a tester emit a certificate which any third party can verify cheaply, an order-finding mode (`-order`), a Fermat/GF/xGF divisibility search (`-divides`), and a batch driver (`-batch`).

Two optional build flavors are gated by `#define`s in `src/version.h`:
- `BOINC` — distributed-computing wrapper
- `NETPRST` — networked client

When defined, each adds a runtime subcommand to `main()`'s option parser (`-boinc` at `prst.cpp:199-201`, `-net` at `prst.cpp:202-204`) which `exit()`s into `boinc_main` / `net_main` (in `boinc.cpp` / `net.cpp`). They reuse the same test machinery; they don't replace `main()`.

## Repo layout

```
prst/                            ← repository root
├── README.md                    ← readme
├── sample.ini                   ← example config
├── docs/                        ← this folder (PRST application docs)
├── src/                         ← PRST-specific code (~10k LoC)
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
└── framework/                   ← the shared library (submodule → patnashev/arithmetic, documented there)
    ├── arithmetic/              ← Giant, GWArithmetic, lucas, edwards, montgomery, poly …
    ├── gwnum/ gmp/ bow/         ← prebuilt GWnum / GMP / BOINC state-serialization glue
    ├── inputnum.* logging.*     ← number parser; Logging / SubLogging / Progress
    ├── task.* file.*            ← Task / InputTask base classes; File checkpoint abstraction
    └── config.* container.* md5.*  ← option-parser DSL; FileContainer (proof packs); hashing
```

Mental model: **`src/` is "what test do we run and how do we orchestrate it?"; `framework/` is "C++ wrapper of bignum libraries + UI and IO primitives every tool reuses."** The framework lives in a separate repo (`patnashev/arithmetic`) and is consumed by both PRST and other Atnashev-family tools — so its docs live with it (see below).

## End-to-end execution flow

A typical invocation `prst "30006!4-1" -d` flows:

1. **`main()` in `src/prst.cpp:43`** — install signal handlers, declare `Options`.
2. **Option parsing** — uses the framework `Config` DSL to parse argv into `Options`, `proof_op`, log level, etc. (`prst.cpp:82` onward.)
3. **Input parsing** — `InputNum::parse` (framework) classifies the candidate as one of `KBNC` / `FACTORIAL` / `PRIMORIAL` / `GENERIC` and stores k, b, n, c, factors.
4. **Trial-division shortcut** — small inputs (`bitlen <= 40`) or `-trial` are handled inline at `prst.cpp:309`.
5. **Progress file wiring** — a framework `File` for `prst_<fingerprint>*.param` is opened so progress survives restarts (`logging.file_progress(...)` at `prst.cpp:352`).
6. **`Run::create()` (`prst.cpp:463`)** — the dispatcher, in code order:
   - `-order` → `Order`
   - `-divides` → `FermatDivisor`
   - `-proof cert` → reuse the `Proof` instance directly
   - forced `-fermat`, generic numbers, or `|c| ≠ 1` → `Fermat::FERMAT`
   - Proth-form `k*2^n+1` with `k < 2^n` → `Fermat::PROTH`
   - call `expand_factors()`; if not half-factored → `Fermat::AUTO` with a "Not enough factors" warning
   - `c == 1` and half-factored → `Pocklington` (when the number of factors is small) else `PocklingtonGeneric`
   - `c == -1` and half-factored → `Morrison` (when the number of factors is small) else `MorrisonGeneric`
7. **Optional proof wrapping** — if `-proof save|build` is requested, the `Run` is wrapped in a `Proof` that delegates the inner Fermat test (`prst.cpp:358-381`).
8. **GWState setup** — initialize GWnum library, handle all kinds of exceptions (like header/lib version mismatch).
9. **`run->run(gwstate, file_checkpoint, file_recoverypoint, logging)`** — the polymorphic test entrypoint at `prst.cpp:438`. The selected test class drives one or more `Task` instances, writes checkpoints to `prst_<fingerprint>.ckpt`/`.rcpt`, and emits result lines via `logging.result(...)`.
10. **Exit code** — `PRST_EXIT_PRIMEFOUND` (2) / `PRST_EXIT_NORMAL` (0) / `PRST_EXIT_FAILURE` (1), defined in `src/prst.h`.

## The Run hierarchy

Introduced in upstream commit `4764aa9 Refactoring, base class for all tests.` All tests now derive from a common `Run` (`src/prst.h:75`):

```
Run                         (abstract base — name, fingerprint, success/prime flags, run())
├── Fermat                  (Fermat / Proth / Pocklington-base / AUTO; uses Exp)
│   └── Pocklington         (small-factor Pocklington; reuses Fermat::run)
├── PocklingtonGeneric      (FactorTree-driven; runs Exp sub-tasks under a SubLogging)
├── Morrison                (small-factor Morrison; uses LucasMul)
├── MorrisonGeneric         (FactorTree-driven; runs LucasMul sub-tasks under a SubLogging)
├── Order                   (multiplicative order; uses Exp)
│   └── FermatDivisor       (-divides: searches F/GF/xGF numbers divisible by the input)
└── Proof                   (runs the verification for -proof cert; uses Exp)
    └── Proof+Fermat        (orchestrates a wrapped Fermat for -proof save/build)
```

Two patterns repeat:
- **Small-factor variants** (`Morrison`, `Pocklington`) work directly on the parent `Logging`. They are picked for smooth inputs with small number of factors, defined as `n > 10` and `factors.size() < 10` (`prst.cpp:512, 522`). The conditions are subject to change.
- **`*Generic` variants** (`MorrisonGeneric`, `PocklingtonGeneric`) build a `FactorTree`, walk it, and run dynamically-created sub-tasks under their own `SubLogging _logging`. This is where sub-task time accounting happens, and it is a subtle area (the inner `SubLogging` progress must be propagated to the parent — see the framework's `logging-and-progress.md`).

`Run::create` (a static factory in `prst.cpp:463`) is the single source of truth for which class handles which form.

## The Task layer

Tests don't do the heavy math themselves — they construct one or more `Task` subclasses from the framework + `src/exp.h` / `src/lucasmul.h`, configure them, and call `task->run()`. The `Task` / `InputTask` base classes live in the framework (`task.{h,cpp}`); their restart/checkpoint lifecycle is documented in the framework's `task-lifecycle.md`. This doc covers the PRST-side concrete tasks and how the tests drive them.

Concrete tasks of interest (all defined in `src/`):
- `CarefulExp`, `MultipointExp`, `BaseExp` — exponentiation variants with strong-error-check (Gerbicz/Gerbicz-Li) options. Defined in `src/exp.h`.
- `LucasVMul`, `LucasVMulFast`, `LucasUVMul` — Lucas-sequence multipliers. Defined in `src/lucasmul.h`.
- `Product` — multiplies giants under error-checked arithmetic; constructed once, then `mul(a, b)` (two giants in one call) or `mul(first, last)` per multiplication. Used in MorrisonGeneric/PocklingtonGeneric and FermatDivisor. Defined in `src/exp.h`.

A test class typically:
1. Constructs its tasks in its constructor (using `add_stage` to register expected work for progress reporting).
2. Resumes from `progress().param()` if a test state file exists.
3. Sets up unique state files for individual tasks.
4. Loops: call `task->run()`, inspect `task->result()`, advance the algorithm.
5. Writes the result line via `logging.result(...)` and `logging.result_save(...)`.

## The framework (shared library)

PRST builds on the **arithmetic framework** (`patnashev/arithmetic`, the `framework/` submodule). These
primitives are documented in that repo's `docs/` folder — only their PRST-facing role is summarized here.

### `Logging` — user-facing accounting

Message output, results, factor/log files, progress checkpoints, and staged progress/time accounting.
Tests emit results through `logging.result(...)`; `*Generic` tests run inner work under a `SubLogging`.

(Framework: `logging-and-progress.md`.)

### `InputNum` — the parsed candidate

Stores k/b/n/c, form type, factor lists and pretty-printed text. Most of `Run::create`'s branching keys off this.

Factor lists may contain factorials and primorials to speed up parsing. If a test needs to iterate through
individual factors, `expand_factors()` is called first.

(Framework: `inputnum-parsing.md`.)

### `File` — checkpoint and progress storage

PRST (`File::FILE_APPID = 4`) can write several temporary files on disk for each test. Files have generated names with `prst_` prefix, although collisions are possible.

The progress file (`*.param`) is a text `key=value` file. Can be used to obtain `time_total` by an external observer.

Checkpoints are saved in `*.ckpt` and optionally `*.rcpt` files (recovery point of strong checks). Checkpoints can have children (`.1`, `.2`, ...) for sub-tasks.

`LLR2File` is a variant with a different on-disk format used for some compatibility cases (`prst.cpp:336`).

Fingerprint prevents restoring state from a different test run. Therefore, it MUST be changed whenever a parameter which affects calculations changes. Such parameters are:
- input number.
- type of test.
- base of a Fermat test.
- Lucas chain parameters.
- number of proof files.

(Framework: `state-serialization.md`.)

### `FileContainer` — the `.pack` bundle

`-proof … pack` option writes all files necessary for `-proof build` into a single container file. Used only for tests with `-proof` option.

(Framework: `container-format.md`.)

### `Config` — the option-parser DSL

Turns `argv` or `*.ini` into `Options`.

(Framework: `config-dsl.md`.)

### `Options` — user-facing knobs

Bag of `std::optional`s populated by the `Config` DSL. Test classes consult fields like `Check`, `CheckStrong`, `StrongCount`, `FermatBase`, `AllFactors`, `ProofPointFilename`.

Also contains fields which configure GWState: thread count, spin threads, instruction set (SSE2/AVX/FMA3/AVX512F), FFT size hints, safety margin. `Run::create` may also change these fields as a side-effect.

### `GWState` — the math runtime

Owned by `main()`. Can't be intialized before `Run::create` sets the required parameters in the `options`, which are propagated by `options.configure(gwstate)`. `input.setup(gwstate)` intializes GWState, picks the actual FFT and reports it via `gwstate.fft_description`. `gwstate` can be reused after a call to `gwstate.done()`.

(Framework: `arithmetic-foundation.md`.)

## Build

- Windows: `Build commands.txt` invokes MSBuild on `src/win64/PRST.sln` (Release/x64, Win SDK 10.0.26100.0). `BOINC`/`NETPRST` builds are separate solution configs.
- Linux: `src/linux64/Makefile` (and `Makefile.boinc`).
- macOS: `src/mac64/Makefile`.

GWnum + GMP are prebuilt static libs shipped in the framework submodule (`framework/gwnum/{linux64,mac64,win64}`, `framework/gmp/win64`). Bumping GWnum means bumping the framework submodule pointer.

### Version

PRST is a scientific software, and an important aspect of any science is reproducibility. That's why all serious researchers log versions of all software used. If a bug is found in a software, only affected results can be redone. That's why it's important to have a strict binding of reported version to code in the repository. The best solution is to have automatic version increment before each build. All forks of PRST should also follow this rule. They should indicate they're not PRST (for example, by a different major/minor version in a range not overlapping with PRST or other forks). And one should never publish two binaries with the same version number, even if it's a quick build for a friend to perform a single test. Such builds have a tendency to stick for years.

## Working with this codebase — practical notes

1. **Read the framework before touching it.** The framework is a shared library used by other tools; its surface area is tighter than `src/`. Adding a method to the framework's `logging.h` is a much higher-cost change than modifying a caller in `src/morrison.cpp` — and lives in a separate repo/PR (`patnashev/arithmetic`).
2. **`Run::create` is the dispatcher.** Any new test mode plugs in there. Any change to which class handles which input form lives there.
3. **Result lines are user-visible contract.** The `logging.result(...)` strings (`is prime!`, `is not prime`, `time: %.1f s.`) are what mersenneforum users / BOINC parsers consume — don't change wording lightly.
4. **Checkpoint format is on-disk contract.** `TaskState::read/write` and the `File` fingerprint scheme define resumability. Bumping any of them needs a version bump and ideally backward-compat reads.
5. **Single-repo fixes preferred.** Coordinated `prst` + `framework` PRs are reviewable but unloved by upstream. Default to "can I express this with the existing framework API?" before reaching for a new one.
6. **Bug reports come from mersenneforum.org.** Search the forum threads for the exact result-line snippet a user is complaining about — that's typically the fastest path to the failing code path.

## Companion deep-dives in this folder (PRST application)

- **`run-hierarchy.md`** — Fermat / Pocklington(Generic) / Morrison(Generic) / Order at the call-graph level; the `Run::create` dispatcher; the small-factor vs. `*Generic` boundary; `FactorTree` construction and walk; the result-line contract.
- **`proof-system.md`** — the `-proof save / build / cert` pipeline; `Proof` as both `Run` and Fermat-wrapper; `ProofSave` / `ProofBuild` tasks; the binary product-tree schedule; roots-of-unity and security-seed defenses; the `.pack` container.
- **`abc-batch.md`** — `-batch` as an alternate entry point; `CandidateSource` + `parse_batch_file`; the raw / ABC / ABCD / ABC2 formats and `ABCTemplate`; the `batch_main` driver loop; two-level progress/resume (`cur`/`primes`/`composites`); the `-stop on {error,prime,composites,kprime}` conditions.
- **`boinc-and-net.md`** — the `BOINC` / `NETPRST` build flavors; `boinc_main` / `net_main` as `exit()`-into-sibling entry points; `BoincLogging` (control hooks: `state_save_flag`/`heartbeat`/trickles) vs. `NetLogging` (data hooks); `NetFile`/`LLR2NetFile` checkpoint shipping over HTTP; the reduced inline dispatcher; **`net.cpp` is stale vs. the run() refactor**.
- **`test-harness.md`** — `PRST -test`: the `Test`/`DeterministicTest` classes and the `test.data` vector arrays (`res64`/`cert64` oracle); `Test::run`'s SAVE→BUILD→CERT cross-check (incl. the cross-FFT consistency check); the subsets; `RootsTest` (proof-forgery defenses) and `ABCParserTest`.
- **`exponentiation-algorithms.md`** — the exp/Lucas `InputTask`s that do the modular work: `BaseExp`/`MultipointExp` (point schedule, sliding window, smooth-vs-non-smooth) and `StrongCheckMultipointExp`/`LucasUVMul` (the Gerbicz / Gerbicz-Li block check — `L≈√iters`, the `D` check-accumulator, recovery-point rollback); how `Fermat::run` picks the task; the `State` classes.
- **`checkpoints.md`** — what PRST persists: the `TYPE` registry (including the TYPE 5 strong-check placeholder), the state class tree, the `.ckpt`/`.rcpt` strong-check handshake, and the LLR2 on-disk compatibility munging. The generic wire format lives in the framework's `state-serialization.md`.
- **`math-and-theorems.md`** — the *why*: Fermat PRP vs. the deterministic Proth (iff) and Pocklington/BLS (`N−1`, `F>√N`) and Morrison/BLS-Theorem-14 (`N+1`, Lucas) criteria; the half-factored threshold; Gerbicz / Gerbicz-Li error checking; IBDWT multiplication; the theorem→code map. A reading guide tying each test's correctness to its code and to the authoritative references (BLS 1975, eprint 2023/195, the `mult_*.pdf`s) — proofs cited, not reproduced.

### Framework deep-dives (in the `patnashev/arithmetic` repo, `docs/`)

The shared-library subsystems PRST consumes are documented alongside their code: `task-lifecycle.md`,
`logging-and-progress.md`, `inputnum-parsing.md`, `state-serialization.md`, `arithmetic-foundation.md`,
`curves-and-polynomials.md`, `config-dsl.md`, `container-format.md`. Start from that repo's
`lay-of-the-land.md`.

## Coverage status

**Every in-scope PRST subsystem has a deep-dive** (the seven above); the framework subsystems are covered
by the companion series in the `patnashev/arithmetic` repo. A durable caveat worth keeping in mind: **a
"complete" claim has to survive a grep one level deeper, and lower-traffic, lower-level code is *more*
worth documenting, not less.**

What's *intentionally* left without its own doc on the PRST side:

- **`main()`'s top half** — option parsing, the `bitlen ≤ 40` trial-division shortcut, and progress-file wiring before `Run::create` — is walked step-by-step in "End-to-end execution flow" above; that counts as covered.
- **The full mathematical proofs** (BLS Theorem 14, the Gerbicz-Li soundness bound, the IBDWT error bound) live in the cited references; `math-and-theorems.md` is the reading guide to them, not a re-derivation.

A caution for maintainers: the docs reflect the code at the time of writing and several embed verified line numbers and "this is the only live caller" claims — re-grep before trusting any such claim. In particular, citations into `framework/...` files point at code that lives in the `patnashev/arithmetic` repo.

## Cross-cutting open questions (parked)

- **Thread-safety boundary.** PRST's threads live inside GWnum's FFT layer; everything above is single-`Task` on the main thread. If multi-`Task` parallelism is ever added (e.g. parallel inner factors in `MorrisonGeneric`), the entire framework `Logging`/`Progress` subsystem needs locking. See the framework's `logging-and-progress.md` and `task-lifecycle.md`. Each thread needs its own `gwstate` instance obtained via `gwstate.clone()`. This is implemented in `patnashev/prefactor` utility.
- **`_smooth` exponentiation path.** The smooth (`b^n` by repeated squaring / windowed powering) vs. non-smooth (sliding-window over the full exponent) split is documented in `exponentiation-algorithms.md` §1, §2.

## How to write the next deep-dive

The template is the framework's `logging-and-progress.md` / `task-lifecycle.md`. Sections, in this order:

1. One-paragraph intro + source file list + cross-references.
2. Class hierarchy / data-model walkthrough.
3. Annotated source for the central function/method (verbatim from the codebase, with inline comments).
4. Field/method reference, every line annotated — no bare lines.
5. Lifecycle / sequence section.
6. Common patterns.
7. Pitfalls, anchored to real bugs where possible.
8. Quick-reference table.
9. Open questions / explicit non-coverage, with pointers to the relevant other doc(s).

When a new deep-dive lands, add its bullet to "Companion deep-dives in this folder" above (or the
framework repo's index if it documents a framework subsystem), and remove anything it now covers from a
sibling doc's "open questions."
