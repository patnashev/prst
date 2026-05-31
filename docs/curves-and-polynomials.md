# Curves & polynomials — deep dive

The specialized arithmetic that sits *beside* the bignum substrate: elliptic-curve point groups (Edwards, Montgomery), Lucas sequences, and polynomial multiplication — all built on top of `GWNum` (their coordinates *are* `GWNum`s). This is general-purpose framework arithmetic (the `patnashev/arithmetic` submodule serves other tools too), so **most of it is dormant in PRST's own live path** — a fact worth stating plainly up front (verified by grepping callers):

- **`LucasV` / `LucasUV` — LIVE.** The Morrison test's `LucasVMul`/`LucasUVMul` tasks (`exponentiation-algorithms.md`) drive this arithmetic (`lucasmul.cpp`, `prst.cpp`), and its small-prime chain indices come from this layer's precomputed differential-addition-chain table. This is the one part PRST actually runs.
- **Edwards / Montgomery curves — DORMANT.** The ECM cofactor-factoring code is present but its only PRST caller — the probe in `factorize` — is **commented out** (`inputnum.cpp:157-268`; see `inputnum-parsing.md` §6/§8). No live code path instantiates `EdwardsArithmetic`; today it's exercised only by the standalone `arithmetic/test.cpp`.
- **`PolyMult` / `Poly` — DORMANT (no caller anywhere).** The proof product tree does **not** use it — `ProofBuild` multiplies its residues as `Proof::Product` *giants* via `gw` (`proof.cpp:552,663`), not the polynomial engine. `PolyMult` has **no caller anywhere in the repo** — its only mention outside `poly.{h,cpp}` is the `friend class PolyMult;` declaration (`arithmetic.h:229`), not a call — and `poly.cpp` is not even compiled into either build (absent from the linux64 `Makefile` object list and the win64 `PRST.vcxproj`).

> **Scope and honesty.** This is the largest subsystem by raw volume — `group.cpp` alone is ~42k hand-written lines, and the `.cpp` bodies for the curve/poly/Lucas formulas were **not** read line-by-line for this doc. What follows is an **interface + usage map**: the class hierarchy, which arithmetic each consumer uses, and how chains are built — grounded in the headers. The actual addition/doubling/ladder *formulas* and the polynomial-FFT internals are the mathematics, deferred to `math-and-theorems.md` and the framework's `mult_*.pdf`s. Treat any formula-level claim as out of scope here, not as verified.

Source files: `framework/arithmetic/group.{h,cpp}` (the group framework + chain helpers), `edwards.{h,cpp}`, `montgomery.{h,cpp}`, `lucas.{h,cpp}`, `poly.{h,cpp}`.

Prereqs / companions: `arithmetic-foundation.md` (every coordinate is a `GWNum`; these follow the same `FieldArithmetic`/`FieldElement` strategy pattern, as `GroupArithmetic`/`GroupElement`), `exponentiation-algorithms.md` (`LucasVMul`/`LucasUVMul` are the *tasks*; `LucasV`/`LucasUV` here are the *arithmetic*), `inputnum-parsing.md` (the ECM factoring consumer), `proof-system.md` (the product-tree consumer), `math-and-theorems.md` (the formulas).

## 1. Two kinds of group arithmetic

Everything here specializes one of two abstract bases in `group.h`, mirroring the `FieldArithmetic`/`FieldElement` pattern from the bignum layer:

- **`GroupArithmetic<Element>`** (`group.h:13`) — a *full* group: `add`, `sub`, `neg`, `dbl`, and scalar `mul` via a windowed NAF chain (`mul(a, W, naf_w, res)`). You can add any two elements.
- **`DifferentialGroupArithmetic<Element>`** (`group.h:197`) — a *differential* group: `add(a, b, a_minus_b, res)` (addition needs the difference `a−b` as a third input — the Montgomery/Lucas "ladder" property), `dbl`, and scalar `mul` via a differential addition chain. No free `add`; you climb a ladder.

This split is the organizing principle. The five concrete arithmetics:

| Arithmetic | Base | Element | Used for (PRST liveness) |
|---|---|---|---|
| `EdwardsArithmetic` | `GroupArithmetic` | `EdPoint` (X,Y,Z,T) | ECM, twisted-Edwards full-group form — **dormant** (§4) |
| `MontgomeryArithmetic` | `DifferentialGroupArithmetic` | `EdY` (Y,Z) | ECM Montgomery-ladder (y-only) form — **dormant** (test-only) |
| `LucasVArithmetic` | `DifferentialGroupArithmetic` | `LucasV` (V, parity) | **live**: Morrison's V-only chain (`LucasVMulFast`) |
| `LucasUVArithmetic` | `GroupArithmetic` | `LucasUV` (U,V, parity) | **live**: Morrison's Gerbicz-checked chain (`LucasUVMul`) |
| `PolyMult` / `Poly` | — | `Poly` (vector of `gwnum`) | polynomial multiply engine — **no caller anywhere; `poly.cpp` not built** (§4) |

The element types hold their coordinates as `unique_ptr<GWNum>` — e.g. `EdPoint::{X,Y,Z,T}` (`edwards.h:127-130`), `LucasUV::{U,V}` (`lucas.h:264-265`) — so a "point operation" is a handful of `GWArithmetic` mul/add/sub on those `GWNum`s.

## 2. Chain construction — NAF windows and DAC

Scalar multiplication (`point × scalar`) is where the two group kinds diverge, and it's the part this layer most clearly owns:

- **Full groups → windowed NAF.** `get_NAF_W(W, a, res, compress)` (`group.h:10`) decomposes the scalar `a` into a width-`W` non-adjacent form — a signed-digit representation that minimizes nonzero digits — into `res` (a `vector<int16_t>`). `GroupArithmetic::mul(a, W, naf_w, res)` then walks those digits with a precomputed odd-multiples table (the group analogue of `MultipointExp`'s sliding window). Edwards ECM and `LucasUV` use this.
- **Differential groups → DAC.** A differential addition chain (where every `add` has its difference available) is built from `get_DAC_S_d(e, start, end, maxlen)` plus a baked-in table `precomputed_DAC_S_d[]` of length `precomputed_DAC_S_d_len` (`group.h:192-194`). For a small prime `e`, the chain is looked up in the table; for larger ones it's searched. **This is the source of Morrison's `_dac_index`** — `MorrisonGeneric` precomputes a DAC index per small prime factor (`morrison.cpp`) by indexing into `precomputed_DAC_S_d`, and `LucasV`/Montgomery climb the corresponding ladder.

The `DifferentialGroupArithmetic::mul(a, prime, index, res)` overload (`group.h:334`) keys off exactly this: `index < precomputed_DAC_S_d_len` → use the table entry, else compute one via `get_DAC_S_d` (`group.h:347-350`).

## 3. The concrete arithmetics

**`EdwardsArithmetic` / `EdPoint`** (`edwards.h`) — twisted Edwards curves for ECM. `EdPoint` carries extended projective coordinates `(X,Y,Z,T)`; `extend()` materializes `T = X·Y` and projects, `normalize()` reduces to affine. Beyond the group ops it offers the ECM-specific surface: `gen_curve(seed, ed_d)` (generate a pseudo-random curve from a seed), `from_small(...)` (a curve from small integer parameters), `on_curve(a, ed_d)` (membership test), `jinvariant`, `d_ratio`. `ED_PROJECTIVE`/`EDADD_NEGATIVE`/`EDDBL_FOR_EXT_NORM_ADD` are `option` flags into `add`/`dbl`. (`arithmetic/test.cpp` exercises this path: `gen_curve`, `on_curve`, `dbl`/`add`, `mul` with a NAF.)

**`MontgomeryArithmetic` / `EdY`** (`montgomery.h`) — the differential (y-coordinate-only) form, constructed *from* an `EdPoint` (`init(const EdPoint&, EdY&)`) and carrying `(Y,Z)` plus cached `ZpY`/`ZmY`. This is the cheap ladder used for the bulk of ECM scalar multiplication; `add` takes the `a_minus_b` difference, `dbl`, `optimize`.

**`LucasVArithmetic` / `LucasV`** (`lucas.h:9-105`) — the V-only Lucas sequence (differential), parameterized by `negativeQ`. `LucasV` holds `V` (a `GWNum`) and a `parity` bit. This is the arithmetic behind `LucasVMulFast` (the non-checked Morrison chain).

**`LucasUVArithmetic` / `LucasUV`** (`lucas.h:109-268`) — the full `(U,V)` Lucas pair (a full group, so it has `add`/`sub`/`neg` and NAF `mul`), behind the Gerbicz-checked `LucasUVMul`. Notable: a **small-value optimization** — the constructor precomputes a table `_UV_small` of `(U,V,DU)` triples for a small integer `P` (`lucas.h:122-149`), so chains over small `P` avoid full `GWNum` work; `dbl_add_small`/`init_small`/`max_small` use it. Falls back to a `GWNum` `_P`/`_D` for large `P`.

**`PolyMult` / `Poly`** (`poly.h`) — polynomial multiplication over `GWNum` coefficients, wrapping the GWnum `polymult` library (`pmhandle`). `Poly` is a vector of `gwnum` with a `monic` flag and an optional preprocessed `_cache`. `PolyMult` offers a full toolkit: `mul`, `mul_split`, `mul_range`, `mul_twohalf`, `fma_range`, `reciprocal`, `preprocess`(+`preprocess_and_mul`), `eval`. **It has no caller anywhere in the repo** — the proof BUILD product tree multiplies *giants* (`Proof::Product`, `proof.cpp`), not polynomials; nothing in the tree names `PolyMult` except the `friend class PolyMult;` declaration (`arithmetic.h:229`), and `poly.cpp` is not even built (it is absent from the linux64 `Makefile` object list and the win64 `PRST.vcxproj`). It is dormant framework code, available for the framework's other tools.

## 4. Usage: who actually calls this

Only one consumer is live in PRST; the rest is dormant-but-present (grep-verified above). Stated honestly:

1. **Morrison — the one live consumer** (`morrison.cpp` + `lucasmul.cpp`, reached from `prst.cpp`). `LucasVMulFast`/`LucasUVMul` advance `LucasV`/`LucasUV` state; small odd primes are multiplied via the DAC ladder, with `MorrisonGeneric` precomputing the chain index per small prime (the `_dac_index` array) from `precomputed_DAC_S_d`. The Gerbicz-checked variant uses the full `LucasUV` group. This is the part of the layer PRST runs in production.
2. **ECM factoring — dormant.** The elliptic-curve method *would* pick a pseudo-random Edwards curve (`gen_curve`), use the Montgomery ladder for the bulk scalar-multiply, and read a factor from a non-trivial `gcd` of a coordinate with `N` — but the code that does this is the **commented-out** block in `InputNum::factorize` (`inputnum.cpp:157-268`). The live `factorize` is trial-division only (`inputnum-parsing.md` §6). So the curve arithmetic is present and unit-tested (`arithmetic/test.cpp`) but not wired into any live PRST run.
3. **`PolyMult` / `Poly` — dormant, no caller anywhere.** The proof product tree does *not* use this layer: `ProofBuild` compresses proof points by multiplying *giants* (`Proof::Product`, `proof.cpp:552,663` via `gw`), not the `PolyMult` polynomial engine — a natural assumption that the grep refutes. In fact `PolyMult` has **no caller anywhere in the repo**: its only mention outside `poly.{h,cpp}` is the `friend class PolyMult;` declaration (`arithmetic.h:229`), and `poly.cpp` is not even compiled into either build (absent from the linux64 `Makefile` object list and the win64 `PRST.vcxproj`).

Where it *is* used (Morrison, and the dormant curve/poly code), the *coordinates/coefficients* are `GWNum`s in the same `GWState` FFT the test runs in — so this layer is "structured arithmetic on top of the one modular-multiply engine," not a separate number system.

## 5. Common patterns

- **Coordinates are `GWNum`s, ops are `GWArithmetic` calls.** Every `EdwardsArithmetic::add`/`dbl`, every `LucasUVArithmetic::mul` is, underneath, a sequence of `gw().mul`/`square`/`add` on the element's coordinate `GWNum`s — so the round-off / reliability machinery and the FFT config all come from `arithmetic-foundation.md` unchanged.
- **Full vs differential is a real fork, not a style.** A `GroupArithmetic` lets you `add` anything (Edwards, `LucasUV`); a `DifferentialGroupArithmetic` only lets you `add(a,b,a−b)` and so forces a ladder (Montgomery, `LucasV`). Choosing the differential form trades generality for ~half the work — which is why Morrison's cheap path uses `LucasV` and ECM's bulk uses the Montgomery ladder.
- **NAF for full, DAC for differential.** Two distinct scalar-multiply chain builders (`get_NAF_W` / `get_DAC_S_d`+table), picked by which group kind you're in.
- **Small-integer fast paths.** `LucasUVArithmetic`'s `_UV_small` table and the `precomputed_DAC_S_d` table both hard-code chains for small operands to skip full `GWNum` arithmetic.

## 6. Pitfalls

- **This doc maps the surface, not the formulas.** The curve addition/doubling, the Montgomery ladder, the Lucas recurrences, the NAF/DAC construction, and the polynomial FFT all live in `.cpp` bodies (notably ~42k lines of `group.cpp`) not walked here. Don't treat the descriptions above as a substitute for reading the relevant `.cpp` before modifying a formula — and cross-check the math against `math-and-theorems.md`/the PDFs.
- **`EdPoint` coordinates may be null / unextended.** `X/Y/Z/T` are `unique_ptr<GWNum>`; `Z`/`T` can be absent until `extend()`/`normalize()`. Code reading `*T` without ensuring extension will dereference null.
- **Differential `add` needs the *right* difference.** `add(a, b, a_minus_b, res)` is only correct when `a_minus_b` truly equals `a−b` on the ladder; feeding the wrong difference silently computes garbage (no error — it's not checked).
- **`MorrisonGeneric`'s `_dac_index` encodes table-hit vs. computed by *sign*.** For a small prime found in `precomputed_DAC_S_d` it stores the table *position* (positive); for a prime beyond the table it stores `-get_DAC_S_d(...)` — the **negative** of a freshly-computed chain parameter (`morrison.cpp`). It is *not* an "index ≥ len means compute" scheme; the sign is the discriminator (`group.h`'s differential `mul` reads the table for a positive in-range index, else takes the value inline). A wrong value builds the wrong differential chain — a correctness bug, not just a perf one.
- **`Poly` ownership is subtle.** A preprocessed `Poly` holds a `_cache` instead of `_poly`; `size()`/`data()` switch on `_cache != nullptr`. Mixing preprocessed and raw polynomials, or freeing across a `PolyMult`, needs care (the move/copy go through `PolyMult`).

## 7. Quick reference

| You want… | Class / call |
|---|---|
| A fast full elliptic-curve group (ECM) | `EdwardsArithmetic` / `EdPoint` |
| The cheap EC ladder (ECM bulk) | `MontgomeryArithmetic` / `EdY` |
| A random ECM curve from a seed | `EdwardsArithmetic::gen_curve(seed, ed_d)` |
| Morrison's V-only chain | `LucasVArithmetic` / `LucasV` |
| Morrison's checked U-V chain | `LucasUVArithmetic` / `LucasUV` |
| Multiply a point by a scalar (full group) | `get_NAF_W` → `mul(a, W, naf_w, res)` |
| …by a small prime (differential) | `precomputed_DAC_S_d[index]` → ladder `mul` |
| Polynomial multiply (framework engine; no caller anywhere, `poly.cpp` not built) | `PolyMult` / `Poly` |

(The proof product tree multiplies *giants* (`Proof::Product`), not polynomials — see §4. ECM/curve rows above are dormant in PRST today — §4.)

Two group kinds: **full** (`GroupArithmetic`: `add`/`sub`/`neg`/`dbl` + NAF `mul`) vs **differential** (`DifferentialGroupArithmetic`: `add(a,b,a−b)`/`dbl` + DAC `mul`). All coordinates are `GWNum`s.

## 8. Open questions / non-coverage

- **All the formulas.** Edwards/Montgomery point addition & doubling, the j-invariant and `gen_curve` construction, the Lucas `V`/`U` recurrences, the NAF and DAC chain algorithms, and the polynomial-FFT (`polymult`) internals — the actual mathematics — are deferred to `math-and-theorems.md` and `framework/docs/mult_*.pdf`. This doc deliberately stops at the interface.
- **`group.cpp` (~42k lines) is unread here.** A line-level deep-dive of the group arithmetic would be its own (large) effort; what's documented is the public surface from `group.h` plus the consumers. If a curve/ladder bug surfaces, that file is the next stop.
- **ECM strategy and parameters.** *Which* curves/bounds the factorizer tries, stage-1/stage-2 structure, and the (currently commented-out) deeper ECM probe in `inputnum.cpp`'s `factorize` are factoring-strategy questions touched in `inputnum-parsing.md` §6/§8, not re-derived here.
- **`polymult` library boundary.** `PolyMult` wraps GWnum's `polymult` (`pmhandle`, `polymult.h`) — like `gwnum` itself, a prebuilt component understood at its API, not dissected (cf. `arithmetic-foundation.md`'s external-library note).
- **The proof product tree is not in this layer.** `ProofBuild` multiplies `Proof::Product` *giants* via `gw` (`proof.cpp`), not `PolyMult` — its schedule is `proof-system.md` territory and bottoms out in the bignum layer, not here. (Noted because "polynomial product tree" is the natural-but-wrong assumption.)
