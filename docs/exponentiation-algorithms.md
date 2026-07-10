# Exponentiation algorithms — deep dive

This is the missing middle: between the `Task` restart/checkpoint loop (`task-lifecycle.md`) and the number theory of *why* the tests are valid (`math-and-theorems.md`) sit the concrete tasks that actually grind out the modular exponentiation. When `Fermat::run` calls `task->run()`, *this* is the `task` — millions of GWnum squarings, scheduled to hit proof points, optionally wrapped in a Gerbicz / Gerbicz-Li error check. Two parallel families:

- **`exp.{h,cpp}`** — `GWNum` exponentiation, for the Fermat/Proth/Pocklington side. The base hierarchy is `BaseExp` → `CarefulExp`, `MultipointExp` (and `Smooth`/`Fast`/`SlidingWindow`), `StrongCheckMultipointExp` (and `GerbiczCheck`/`LiCheck`/`FastLiCheck`).
- **`lucasmul.{h,cpp}`** — Lucas V/UV-sequence chains, for the Morrison side: `LucasVMul` → `LucasVMulFast`, `LucasUVMul` → `LucasUVMulFast`.

The crux of both is the **strong check** (`StrongCheckMultipointExp`, `LucasUVMul`): Gerbicz's identity (Gerbicz-Li for non-smooth exponents) lets a long exponentiation verify itself periodically and roll back to a recovery point on a detected error — the engineering behind "checked" primality tests.

This doc covers the *task orchestration* — point scheduling, the sliding-window engine, the Gerbicz block structure, the smooth-vs-non-smooth split. It does **not** re-derive *why* the Gerbicz identity holds (→ `math-and-theorems.md`) or how the `LucasUV` arithmetic recurrences work (→ `curves-and-polynomials.md`).

Source files: `src/exp.h`, `src/exp.cpp`, `src/lucasmul.h`, `src/lucasmul.cpp`.

Prereqs / companions: `task-lifecycle.md` (these are `InputTask`s; `execute()` runs inside the restart loop, `commit_execute` writes checkpoints), `arithmetic-foundation.md` (every squaring is a `GWArithmetic::square`/`mul`; intermediate state is a `SerializedGWNum`, final is a `Giant`; `ReliableGWArithmetic` is the round-off layer), `run-hierarchy.md` (which test constructs which task — the `Fermat` ctor's `FastExp`/`GerbiczCheckExp`/`StrongCheck…` selection), `proof-system.md` (the `on_point` callback writes proof points; `MultipointExp::Point` is the proof schedule), `checkpoints.md` (the TYPE registry — 1/2/5/8 here, 9/10/11 in lucasmul — and the `.ckpt`/`.rcpt` file handshake; the wire format is the framework's `state-serialization.md`).

## 1. Class model

**State** (`exp.h:16-57`) — what a checkpoint holds. `BaseExp::State` is abstract (`set`/`to_GWNum`); two concretes:
- `StateValue` (TYPE 1) — an exact `Giant`. Used for the **final** result and for values that must survive (`commit_execute` writes this at `iteration == iterations()`).
- `StateSerialized` (TYPE 8) — a `SerializedGWNum`. Used for **intermediate** checkpoints (cheaper, FFT-domain). `commit_execute` picks `StateValue` at the last iteration, `StateSerialized` otherwise (`exp.h:67-74`).

**Smooth vs. non-smooth** (`_smooth`, `exp.h:78-83`). The single most important fork:
- **smooth** — the exponent *is* `b^n` (a base raised to a power; "smooth" because its factorization is just `b`'s). Computed by repeated squaring (`b==2`: one square per step) or windowed powering (`b>2`: `sliding_window(b^L)` per block). `b()` returns the base.
- **non-smooth** — an arbitrary exponent `_exp` (e.g. `N-1`), computed by sliding-window over its bits, starting from `X0`/`x0`, optionally folding the base in via `setmulbyconst`.

**The exp task tree** (`exp.h`):

| Class | Mode | Role |
|---|---|---|
| `CarefulExp` | exact | short/exact legs via `gw.carefully()`; `cost = bitlen·1.5` |
| `MultipointExp` | either | the workhorse: exponentiate, stopping at each `Point` to fire `on_point` |
| `SmoothExp` / `FastExp` / `SlidingWindowExp` | smooth / non-smooth(`x0`) / non-smooth(`X0`) | single-point convenience wrappers |
| `StrongCheckMultipointExp` | either | adds the Gerbicz/Gerbicz-Li block check (§5) |
| `GerbiczCheckExp` / `LiCheckExp` / `FastLiCheckExp` | smooth / non-smooth / non-smooth(`x0`) | the checked leaves, with point schedules pre-built in the ctor |

**`MultipointExp::Point`** (`exp.h:141-148`): `{int pos; bool check; bool value;}` — `pos` = iteration to stop at, `check` = a Gerbicz check boundary, `value` = materialize as a `StateValue` (Giant) rather than serialized. The `_points` vector is the proof-point / check schedule.

**The Lucas task tree** (`lucasmul.h`) mirrors it: `LucasVMulFast` (plain V-chain, `State` TYPE 9) is the analogue of `MultipointExp`; `LucasUVMul`/`LucasUVMulFast` (`State` TYPE 10, `StrongCheckState` TYPE 11) is the Gerbicz-checked analogue of `StrongCheckMultipointExp`, working over `LucasUV` arithmetic instead of `GWNum` squaring.

Two more `InputTask`s live in these headers outside the exponentiation trees: `Product` (§6) and the proof tasks `ProofSave`/`ProofBuild` (`proof-system.md`).

## 2. `MultipointExp::execute` — the central method

`exp.cpp:182-272`. The point loop: resume from the checkpoint (or seed `X0`/`x0`), then for each scheduled `Point` advance `X` to `pos` by one of three strategies, fire the callback, checkpoint. Annotated:

```cpp
void MultipointExp::execute()
{
    if (result() != nullptr) return;                      // already done
    len = _exp.bitlen() - 1;
    if (_x0 > 0) gw().setmulbyconst(_x0);                 // fold the base into each square
    if (state() == nullptr && !smooth()) { /* seed X = X0 / x0, checkpoint iteration 0 */ }
    else { i = state()->iteration(); state()->to_GWNum(X()); }   // resume
    if (i < 30) gwset_carefully_count(gw().gwdata(), 30 - i);     // first 30 ops careful (FFT warm-up)

    for (next_point = …; next_point < _points.size(); next_point++)
    {
        if ((smooth() && b() == 2) || (!smooth() && _x0 > 0))     // strategy A: bit-by-bit squaring
            for (; i < _points[next_point].pos; i++) {
                (i < iterations()-30 ? gw() : gw().carefully())
                    .square(X(), X(), (!smooth() && _exp.bit(len-i-1) ? GWMUL_MULBYCONST : 0)
                                      | GWMUL_STARTNEXTFFT_IF(…));   // mulbyconst injects an exp bit
                if (i+1 != pos) commit_execute(i+1, X());            // serialized checkpoint
            }
        else if (smooth())                                          // strategy B: windowed powering
            { tmp = b(); tmp.power(pos - i); sliding_window(tmp); i = pos; }
        else if (!_X0.empty())                                      // strategy C: windowed over _exp bits
            { slide(_exp, i, pos, true); i = pos; }

        check();                                                    // ReliableGWArithmetic round-off gate
        tmp_state = State::cast(_points[next_point].value, _tmp_state);  // Giant or serialized
        tmp_state->set(i, X());
        if (_on_point != nullptr) { … _on_point(next_point, tmp_state, *_logging); }   // proof point
        _tmp_state.swap(_state); on_state();                        // checkpoint at the point
    }
    if (i < iterations()) { /* final tail multiply: X *= _tail */ }
    done();
}
```

The three strategies: **(A)** plain binary exponentiation — for smooth `b==2` (just square `n` times) and non-smooth-with-`x0` (where each set exponent bit is injected by `GWMUL_MULBYCONST`); **(B)** smooth `b>2` — raise to `b^(power)` per segment via `sliding_window`; **(C)** non-smooth-with-`X0` — slide over `_exp`'s bits. The first 30 ops run carefully (`gwset_carefully_count`) to let the FFT settle, and the last 30 run carefully too (`i < iterations()-30 ? gw() : gw().carefully()`). The trailing `_tail` multiply applies the small correction factor the test framework computed (e.g. for `c ≠ 1`).

**Sliding window** (`slide_init`/`slide`/`sliding_window`, `exp.cpp:274-351`): precompute the odd-power table `_U[0..2^(W-1))`, choosing window width `W` by a cost heuristic (`slide_init`'s loop weighs table cost `2^(W-1)` against per-bit savings `len·(1+1/(W+1))`), capped by `_W` (default 5) and `_max_size`. `slide` then walks the exponent bits, squaring and multiplying by the right table entry per window.

## 3. Field & method reference

The `init_*` family selects the mode and seed (each test's ctor calls exactly one):
- `init_smooth(input, gwstate, file, logging[, tail])` — smooth mode; exponentiate `b^n`.
- `init_small(…, uint32_t x0[, tail])` — non-smooth, integer start `x0` (folded via `setmulbyconst`).
- `init_giant(…, Giant X0[, tail])` — non-smooth, full `Giant` start.
- `StrongCheckMultipointExp` adds a `file_recovery` argument (the recovery-point file) and an optional `tail_inv`.

| Member | Meaning |
|---|---|
| `_points` / `points()` | the stop schedule (`Point{pos,check,value}`); built in the leaf ctors |
| `_on_point` | callback fired at each point (proof-point writer; `proof-system.md`) |
| `_W` (=5), `_max_size` | sliding-window width cap and table-size cap |
| `_tail` / `tail()` | trailing multiplier applied after the last point |
| `_x0` / `_X0` | non-smooth start (small / giant) |
| `smooth()` / `b()` / `exp()` | mode + base + exponent accessors |
| `cost()` | estimated op count. Two consumers: each stage's `cost()` is registered via `progress().add_stage(...)` (`run-hierarchy.md`; e.g. `fermat.cpp:298-299`, `pocklington.cpp:160`), which drives the **progress bar**, and the summed `progress().cost_total()` is printed as the **`complexity = %d`** field in the info header the tests log (`fermat.cpp:326-328`, `pocklington.cpp:60`, `morrison.cpp:219`, `proof.cpp:364`) |
| `result()` | the final `Giant`, or `nullptr` until `iteration == iterations()` |
| `commit_execute(i, …)` | checkpoint: `StateValue` at the end, `StateSerialized` mid-run |

`StrongCheckMultipointExp` adds `_L`/`_L2` (block params), `R()`/`D()` (the recovery and check accumulators), `state()` (the recovery point) vs. `state_check()` (the within-block `StrongCheckState`), and `Gerbicz_params`. `LucasVMulFast` builds its chain with `mul_giant`/`mul_prime`; `LucasUVMul` takes `(exp, count)` and computes `_L`/`_L2` from its own `Gerbicz_params` — which rounds the other way than the exp side's (see §5).

## 4. Lifecycle: which task runs

`Fermat`'s constructor (see `run-hierarchy.md` §4.1) picks the task from the check level and form:
- no strong check → `FastExp` (smooth/`x0`) or plain `MultipointExp`;
- strong check, smooth → `GerbiczCheckExp`; non-smooth → `FastLiCheckExp`/`LiCheckExp`;
- with a proof → `MultipointExp`/`StrongCheckMultipointExp` carrying `Proof::on_point` and the proof-point schedule.

`PocklingtonGeneric`/`MorrisonGeneric` choose per tree node similarly (`CarefulExp` for tiny exponents, `Li`/`SlidingWindow` deeper). Once constructed, `task->run()` enters the `task-lifecycle.md` loop: `setup()` allocates the `GWNum`s, `execute()` runs the point loop, checkpoints land via `commit_execute`/`on_state`, and `result()` yields the final `Giant`. The `Point` schedule does double duty: proof points (where `on_point` writes a `ProofSave` point) and Gerbicz check boundaries (`check == true`).

## 5. The Gerbicz / Gerbicz-Li check

The strong check turns a blind `n`-squaring into a *self-verifying* one. The `n` iterations are divided by the proof points into `count` segments (≈ `n/count` each); `Gerbicz_params` (`exp.cpp:385-400`) is called **per segment** (`iters = n/count`, e.g. `exp.h:444`), choosing `L ≈ √iters` and `L2` = the largest multiple of `L` near `iters` — so `L2 ≈ iters ≈ L²`, and one segment is one Gerbicz block. Two running accumulators over that block:
- `X` — the main value, squared each step;
- `D` — the **check** accumulator: every `L` iterations `D *= X` (`exp.cpp:598-599`), so ≈ `L2/L ≈ L` updates per block.

At the block boundary (`exp.cpp:669-674`) the accumulated `D` is combined with the current `X` and checked against an independent recomputation (square `L` times / `sliding_window(b^L)` / the Li reconstruction) — the Gerbicz identity. A mismatch is caught like any error-check failure and the task **restarts from the recovery point** `_state_recovery` (the last *verified* recovery `State`, written by `write_state`, `exp.cpp:467-475`), redoing the block. That recovery `State` is materialized as a `StateValue` (exact `Giant`) only at value-materialize points (`_points[…].value == true`) and as a `StateSerialized` (FFT-domain) otherwise — its concrete type is chosen by `State::cast`/`read_file` (`exp.cpp:13-67`), so it's a generic `State`, not always a `StateValue`. The within-block position is the `StrongCheckState` (TYPE 2: `recovery` index + serialized `X` + `D`), so an interrupted run resumes mid-block, not just at block boundaries. If `L > 1` and a check point falls inside the remaining distance, `L`/`L2` are halved adaptively (`exp.cpp:557-566`) so the final partial block is still checked.

**Gerbicz vs. Gerbicz-Li.** Smooth exponents (`b^n`) use plain **Gerbicz**; arbitrary exponents (`N-1`, the `LiCheckExp` path) use **Gerbicz-Li**, the Pietrzak/Li generalization for non-smooth exponents — same block machinery, different identity (the `_logging` line prints `Gerbicz-Li check enabled` for non-smooth, `exp.cpp:434`). The `Li` variant's cost includes the extra `log2(L2/L)` term (`exp.cpp:407,420`). `StrongCheckMultipointExp` serves *both* variants because of its lineage: it grew out of the original smooth, Gerbicz-only code, and that legacy is where its edge cases live — the downward `L`/`L2` adaptation at the end of a block (`exp.cpp:552-566`) and the ragged-tail Li reconstruction (`exp.cpp:696-735`).

**The Lucas analogue.** `LucasUVMul::execute` (`lucasmul.cpp`) runs the same block-check structure — `_L`/`_L2` from `Gerbicz_params`, a recovery point (`State`, TYPE 10) and within-block `StrongCheckState` (TYPE 11) — but each "step" advances a `LucasV` pair (`Vn`,`Vn1`) by the chain recurrence (`lucasmul.cpp:278-279`, via `LucasVArithmetic`) rather than squaring a `GWNum`. The `LucasUV` objects are the Gerbicz accumulator/state `X`, `R`, and `D` (`lucasmul.h:201-203`); each step the advanced `LucasV` pair is folded into the `LucasUV` accumulator `D`. The recurrence arithmetic itself belongs to `curves-and-polynomials.md`/`lucas.cpp`; here it's "the same Gerbicz handshake over Lucas state."

**LucasUVMul is the reference strong-check implementation.** Historically it is the *most recent* class implementing the strong check, and it implements **only Gerbicz-Li** (its log line is hardcoded, `lucasmul.cpp:196`) — which lets it be simpler where the exp side is crufty. In particular it has **no end-of-block `L2` adaptation**: its `Gerbicz_params` rounds `L2` *up* to a multiple of `L` (`lucasmul.cpp:178-190`, vs. the exp side's round-*down*, `exp.cpp:385-400`), and the surplus positions are read as **leading zeroes in the exponent** (`substr` past `bitlen()` returns 0 — `giant.cpp:947-958`), so `(_L, _L2)` stay fixed for the whole run and every block is uniform. Performance is the same either way; the newer code is just clearer, with fewer edge cases. **Any new strong-check implementation should look at `LucasUVMul` first**, not `StrongCheckMultipointExp`.

**The first/last-30 careful ops.** Both `execute`s force `gw().carefully()` for the first 30 iterations (FFT warm-up) and the last 30 (so the final residue is exact) — a reliability detail independent of the Gerbicz check.

## 6. `Product` — error-checked giant multiplication

`class Product : public InputTask` (`exp.h:544`). Not an exponentiation: it multiplies giants under the same restartable, error-checked `Task` machinery as everything else here. Constructed once per run scope — `Product Pr(&input, &gwstate, &logging);` — with error checking forced on (`set_error_check(false, true)`) and no checkpoint file (snapshots are in-memory `StateValue`s, enough for the reliable-arithmetic restart path; an aborted product recomputes). Two entry points:

- `Giant mul(iterator first, iterator last)` — folds a range of giants into one product, one operand per iteration, committing via `commit_execute<BaseExp::StateValue>`;
- `Giant mul(Giant& a, Giant& b)` — two giants in one call.

Each `mul()` clears `_values`, resets `_state`, sets `_iterations` to the operand count, calls `run()`, and returns `std::move(result())` — one instance serves many independent products. `execute()` resumes at `state()->iteration()` and multiplies with `gw().carefully().mul`.

Callers: the `*Generic` GCD-batch products (`morrison.cpp:318`, `:528`; `pocklington.cpp:109`, `:376`) and `FermatDivisor`, which uses `mul(a, b)` to assemble composite bases from prime-base values (`run-hierarchy.md` §4.7).

## 7. Pitfalls

- **Smooth vs. non-smooth changes the meaning of `_exp`.** In smooth mode `_exp` is the *base* `b` (and `b()` aliases it); in non-smooth mode `_exp` is the full exponent and `b()` returns a null reference. The `init_*` family and the `GWASSERT(smooth()…)` guards exist to keep these from being mixed; calling `init_giant` on a smooth task (or `b()` on a non-smooth one) is a bug, not a no-op.
- **`commit_execute` chooses the State type by iteration.** Intermediate checkpoints are `StateSerialized` (FFT-domain, cheap); only the final one is a `StateValue` (exact `Giant`). Code reading `result()` mid-run gets `nullptr` until `iteration == iterations()` — that's the signal "not done," not an error.
- **The recovery point and the check state are different files/records.** `state()` is the recovery point (`_state_recovery`, written to `file_recovery`); `state_check()` is the within-block `StrongCheckState` (written to `file`). A Gerbicz rollback restores `X`/`D` from the *recovery* point, discarding the in-block progress. Confusing the two when touching checkpoint logic corrupts the rollback.
- **Window width is a tuned heuristic, capped by `_W`/`_max_size`.** `slide_init` balances table size against per-bit cost; `_W` (default 5) and `_max_size` cap it. Changing these is a performance knob, but `_max_size` also bounds memory (the `_U` table is `2^(W-1)` GWNums).
- **`L`/`L2` adapt downward near the end — on the exp side only.** Don't assume every block is exactly `L2` iterations — the tail block shrinks (`exp.cpp:557-566`) so the last check still happens. A "check never fired" bug usually traces to this adaptation. `LucasUVMul` never adapts: it pads the exponent with leading zeroes instead (§5), so its `_L`/`_L2` are trustworthy constants.
- **`GWMUL_MULBYCONST` is how the base enters a non-smooth squaring.** In strategy (A), the exponent bits aren't applied by separate multiplies — each set bit toggles `GWMUL_MULBYCONST` on that square (`setmulbyconst(_x0)` was set once up front). Miss this and the exponentiation looks like it ignores `_exp`.

## 8. Quick reference

| Want | Class |
|---|---|
| Exact short exponentiation | `CarefulExp` |
| Error-checked product of giants | `Product` (`mul(a, b)` or `mul(first, last)`) |
| Plain `a^(b^n)` (smooth) | `SmoothExp` / `MultipointExp(smooth)` |
| Plain `a^exp` (non-smooth), small base | `FastExp` |
| Plain `a^exp`, giant start | `SlidingWindowExp` |
| Checked smooth exponentiation | `GerbiczCheckExp` |
| Checked non-smooth | `LiCheckExp` / `FastLiCheckExp` |
| With proof points | any `MultipointExp` with `_points` + `on_point` |
| Lucas V-chain (Morrison) | `LucasVMulFast` |
| Checked Lucas U-V chain | `LucasUVMul` / `LucasUVMulFast` |

Gerbicz block (per proof-point segment, `iters = n/count`): `L ≈ √iters`, `L2 ≈ iters ≈ L²`, `D *= X` every `L` steps, verify `D` against a recomputation at the block end, restart from the recovery point on mismatch. State TYPEs: `1` final Giant, `8` intermediate serialized, `2` exp strong-check, `9`/`10`/`11` the Lucas trio.

## 9. Open questions / non-coverage

- **Why the Gerbicz / Gerbicz-Li identity is sound.** This doc documents the *block machinery* (L, L2, D-accumulator, rollback); the proof that the identity detects errors with overwhelming probability — and the Li generalization to non-smooth exponents — is `math-and-theorems.md` (LiCheck → eprint 2023/195).
- **The `LucasUV` recurrence arithmetic.** `lucasmul.cpp` orchestrates the chain; the `LucasUVArithmetic`/`LucasV` operations it calls (`dbl`/`add`/`mul`, the `V`/`U` recurrences) live in `framework/arithmetic/lucas.cpp` → `curves-and-polynomials.md`.
- **The exact Gerbicz_params tuning.** `L = √iters` then a refinement loop (`exp.cpp:394-399`); the `log2b` factor is currently forced to 1 (`exp.cpp:390`, with the original formula commented out). Why that simplification is acceptable is a tuning question, not re-derived here.
- **`cost()` formulas.** The per-class `cost()` estimates (consumed by the progress bar and the logged `complexity` field, §3) encode the window-width and block overheads; they're reproduced as-is, not validated against measured op counts.
- **The smooth exponent's role in test correctness.** *That* the smooth path computes `a^(b^n)` is documented here; *why* that's the right quantity for the Proth/Pocklington/Fermat test (and how `k`, `tail`, `c` fold in) is `run-hierarchy.md` (construction) + `math-and-theorems.md` (justification).
