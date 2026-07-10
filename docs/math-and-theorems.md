# Math & theorems — deep dive

Every other doc answers "what does the code do?" This one answers "*why is the answer correct?*" — the number-theoretic theorems behind each test, and the error-checking and multiplication theory that make a multi-month computation trustworthy. It's the bottom of the stack: the engineering docs (`run-hierarchy.md`, `exponentiation-algorithms.md`, `proof-system.md`) all defer their "why is this valid?" questions to here.

> **What this is — and isn't.** This is a **reading guide and synthesis**, not a from-scratch proof document. The theorems below are standard, well-established results (Proth 1878, Pocklington 1914, Brillhart–Lehmer–Selfridge 1975, Gerbicz 2017, and the Gerbicz-Li checking paper); the *authoritative* proofs live in the cited references, which this doc points to rather than reproduces. The goal is to state each theorem precisely enough to see why PRST's code implements a *correct* test, map it to the implementing code, and tell you where to read the full proof. Formula-level details of the FFT multiplication are in the framework's PDFs and treated as a reference here, not re-derived.

References (all real; the first three are cited directly in the source):
- **BLS** — Brillhart, Lehmer, Selfridge, *New Primality Criteria and Factorizations of 2^m ± 1*, Math. Comp. 29 (1975), `doi:10.1090/S0025-5718-1975-0384673-1` (cited `morrison.cpp:15`). The `N−1` and `N+1` deterministic criteria; "Theorem 14" is the Lucas (`N+1`) one.
- **eprint 2023/195** — `https://eprint.iacr.org/2023/195` (cited `morrison.cpp:14`). The Gerbicz-Li verifiable-computation / error-checking method PRST uses for non-smooth exponents and the proof certificates.
- **`framework/docs/mult_en_20230925.pdf`** (English) / **`mult_ru_20230906.pdf`** (Russian) — the multiplication (IBDWT / FFT) theory.
- Classical: Proth's theorem (1878), Pocklington's theorem (1914), Lucas–Lehmer–Riesel (LLR). Gerbicz's error check (R. Gerbicz, mersenneforum, 2017).

## 1. The probabilistic core — the Fermat PRP test

**Fermat's little theorem.** If `N` is prime and `gcd(a, N) = 1`, then `a^(N−1) ≡ 1 (mod N)`. The converse is *false* — composite `N` can satisfy it (Fermat pseudoprimes; Carmichael numbers pass for all `a`). So `a^(N−1) ≡ 1` proves only that `N` is a **probable prime base `a`** (PRP), never that it's prime.

This is why `Fermat::run` reports "is a probable prime," not "is prime" (`run-hierarchy.md` §4.1), and why PRST prefers a *deterministic* test whenever the factorization of `N∓1` permits one. The Fermat PRP is the fallback (generic numbers, `|c| ≠ 1`, insufficient factors) and the inner engine the deterministic tests build on.

## 2. The `N−1` side — Proth and Pocklington/BLS

These prove primality from a partial factorization of `N−1`.

**Proth's theorem** (for `N = k·2^n + 1`, `k` odd, `k < 2^n`). `N` is prime **iff** there exists `a` with
`a^((N−1)/2) ≡ −1 (mod N)`.
This is a genuine *iff* — a deterministic test, not probabilistic. `Fermat`'s PROTH mode computes that power and checks the residue is `≡ −1` (i.e. `result + 1 ≡ 0` or `N`, `fermat.cpp:417`); `genProthBase` picks a witness `a` with the right quadratic character (`run-hierarchy.md` §4.1).

**Pocklington's theorem / BLS `N−1`.** Write `N − 1 = F · R` with `F` the *fully factored* part. If `F > √N` and, for every prime `q | F`, there is a base `a_q` with
- `a_q^(N−1) ≡ 1 (mod N)`, and
- `gcd(a_q^((N−1)/q) − 1, N) = 1`,

then `N` is prime. (The base may differ per prime `q`; PRST tries successive bases and reuses one across factors when it works, restarting with the next prime base when a factor isn't witnessed — `run-hierarchy.md` §4.2.) (BLS sharpen the bound — `F > N^{1/3}` with auxiliary conditions — but PRST uses the "more than half factored" threshold `F > √N`.) The two `Pocklington` classes implement exactly this: the Fermat congruence first, then the per-prime `gcd` checks, accumulating confirmed factors until `F² > N` (`run-hierarchy.md` §4.2–4.3). **`InputNum::is_half_factored` is this threshold** — it returns true iff the factored product `F` exceeds the cofactor `R` (`inputnum-parsing.md` §5), i.e. `F > √N`.

## 3. The `N+1` side — Morrison and BLS Theorem 14

The dual criterion, using a partial factorization of `N+1` and **Lucas sequences** instead of plain powers. This is what the Morrison test runs (it's ~2× the work of Pocklington, hence the dispatcher prefers `N−1` when both apply).

For a Lucas sequence with parameters `(P, Q)`, discriminant `D = P² − 4Q`, chosen so that `(D/N) = −1`, write `N + 1 = F · R` with `F` factored, `F > √N`. **BLS Theorem 14** gives the primality criterion; the source records the two cases PRST uses (`morrison.cpp:14-18`):
- **`Q = 1`:** `V_l ≡ 2 (mod factor)` for the relevant `l | N+1`.
- **`Q = −1`:** factor 2 is tested for free, and `gcd(U_{(N+1)/q}, N) = gcd(V_{(N+1)/2q}, N)` — the Theorem-14 identity that lets the cheaper `V` sequence stand in for `U`.

`Morrison`/`MorrisonGeneric` pick the smallest valid `P` (via the Kronecker symbol), run the `V`-chain to the required indices, and check these conditions plus the per-factor `gcd`s (`run-hierarchy.md` §4.4–4.5). **LLR** (Lucas–Lehmer–Riesel, for `k·2^n − 1`) is the special case PRST labels "Morrison (LLR) test."

**Lucas sequences themselves.** `V_m` and `U_m` satisfy `V_0=2, V_1=P` and `U_0=0, U_1=1` with the recurrences `V_{m+n} = V_m V_n − Q^n V_{m−n}` (and companions). PRST never materializes the whole sequence — it climbs a *Lucas ladder* (a differential addition chain) to the needed index, which is what the `LucasV`/`LucasUV` arithmetic in `curves-and-polynomials.md` computes.

## 4. Error checking — Gerbicz and Gerbicz-Li

A months-long exponentiation on real hardware *will* hit transient errors (round-off, cosmic-ray bit flips). The strong check makes the computation **self-verifying** cheaply.

**Gerbicz's check** (for smooth exponents — repeated squaring `x_{i+1} = x_i²`, i.e. `a^(b^n)`). A multiplicative invariant is maintained in a check accumulator `D` (update every `L` steps); after a block of `≈ L²` steps, one extra `≈L`-length exponentiation recomputes the invariant, and the two must match. A single check thus validates `L²` squarings at `O(L)` extra cost, so `L ≈ √(block)` minimizes overhead — exactly the `Gerbicz_params` choice (`exponentiation-algorithms.md` §5). A mismatch means an error crept in; the task rolls back to the last *verified* recovery point and redoes the block.

**Gerbicz-Li** (eprint 2023/195) generalizes the identity to **non-smooth exponents** — arbitrary `N−1`, the `LiCheckExp` path — and underpins the proof-certificate soundness. It's the same block machinery (`StrongCheckMultipointExp`, and `LucasUVMul` for the Lucas analogue), with the identity extended so the check still holds when each step folds in an exponent bit rather than a fixed base. The code prints `Gerbicz` for smooth, `Gerbicz-Li` for non-smooth (`exp.cpp:434`, `lucasmul.cpp:196`).

*Why it's sound* — that a passed check implies the block was almost-certainly computed correctly, with a quantifiable error-escape probability — is the content of the references; this doc states the structure, not the probability bound.

## 5. Adjacent theory (pointers)

Two layers have their own homes; the *why* still lives here in spirit:

- **Proof certificates.** The `-proof` system lets a third party re-verify an `O(n)` exponentiation in `O(log n)` work. The soundness is a Fiat–Shamir / sumcheck argument (the Gerbicz-Li lineage, eprint 2023/195) plus a **roots-of-unity** defense against a specific forgery and a security seed. The pipeline and the defenses are in `proof-system.md` §7; `RootsTest` (`test-harness.md` §4) is the executable check that the defenses actually catch forged points.
- **Multiplication (IBDWT).** GWnum multiplies `k·b^n ± c` numbers via Crandall–Fagin's *irrational-base discrete weighted transform*, which folds the modular reduction into the FFT for numbers of that exact form (this is why `GWState::setup` cares so much about the `k`/`b`/`n`/`c` shape and the algebraic `known_factors` trick — `arithmetic-foundation.md` §2). The per-multiply round-off is provably bounded; the `0.4` threshold in `ReliableGWArithmetic` is a conservative cutoff below which the floating-point result rounds to the correct integer. Full theory: `mult_*.pdf`.

## 6. Theorem → code map

| Theorem / method | Condition (informal) | Implemented in |
|---|---|---|
| Fermat PRP | `a^(N−1) ≡ 1` ⇒ probable prime | `Fermat` (`run-hierarchy.md` §4.1) |
| Proth (iff) | `a^((N−1)/2) ≡ −1` ⇔ prime, for `k·2^n+1` | `Fermat` PROTH mode |
| Pocklington / BLS `N−1` | `F>√N` factored + per-prime `gcd` | `Pocklington` / `PocklingtonGeneric` (§4.2–4.3) |
| half-factored threshold | `F > √N` (`F > R`) | `InputNum::is_half_factored` |
| BLS Theorem 14 (`N+1`) | Lucas `V`/`U` conditions, `F>√N` | `Morrison` / `MorrisonGeneric` (§4.4–4.5) |
| LLR | Morrison special case, `k·2^n−1` | `Morrison` "(LLR)" |
| Lucas ladder | `V_m`/`U_m` via differential chain | `LucasV`/`LucasUV` (`curves-and-polynomials.md`) |
| Gerbicz check | `√block` block-verify, smooth | `GerbiczCheckExp` (`exponentiation-algorithms.md` §5) |
| Gerbicz-Li | same, non-smooth + proofs | `LiCheckExp`, `LucasUVMul`, the proof system |
| IBDWT | fast `k·b^n±c` modmul | GWnum (`arithmetic-foundation.md`) |
| multiplicative order | `ord(a) mod N` for fully-factored `N−1` primes | `Order` (`run-hierarchy.md` §4.6) |

## 7. Pitfalls / what's proven vs. assumed

- **PRP ≠ prime.** A Fermat (or Proth-but-not-converse) "probable prime" is not proven; only the deterministic `N∓1` and Proth-iff results are. The result-line wording (`run-hierarchy.md` §7) encodes this distinction deliberately — don't paraphrase "probable prime" as "prime."
- **The deterministic tests assume `F > √N`.** Pocklington/Morrison correctness *requires* the factored part to exceed `√N`; the dispatcher only commits to them when `is_half_factored` holds, falling back to Fermat PRP otherwise (`run-hierarchy.md` §3). Weakening that threshold breaks the proof of correctness, not just performance.
- **Theorem statements here are the standard forms, not PRST's exact parameterization.** PRST chooses specific witnesses (`genProthBase`, the Kronecker-symbol `P` search), specific bounds, and the BLS sharpenings; the precise conditions it checks are in the *code* (the engineering docs) and the *proofs* are in BLS / eprint 2023/195. This doc is the bridge, not the authority — verify against the references before relying on a subtle case (e.g. the exact `Q=±1` Lucas conditions, the `gcd` index arithmetic).
- **Gerbicz soundness is probabilistic in the adversarial sense.** Against random hardware error it's overwhelmingly reliable; the proof-*certificate* soundness against a *malicious* prover is the stronger claim that needs the roots-of-unity defense and the security seed (eprint 2023/195, `proof-system.md` §7) — a passed Gerbicz check alone is not a forgery-proof certificate.
- **The `0.4` round-off threshold is a tuned conservative bound, not a theorem constant.** It's below the provable safe limit with margin; treat it as an engineering parameter justified by the multiplication theory, not a magic number (`arithmetic-foundation.md` §9).

## 8. References & open questions

**Read, in order, for the full picture:**
1. BLS 1975 (`doi:10.1090/S0025-5718-1975-0384673-1`) — the `N−1`/`N+1` criteria and Theorem 14. The foundation for Pocklington and Morrison.
2. eprint 2023/195 — the Gerbicz-Li checking and proof-certificate theory.
3. `framework/docs/mult_*.pdf` — the IBDWT multiplication and round-off theory.
4. A standard reference (Crandall & Pomerance, *Prime Numbers*) for Proth, Pocklington, Lucas sequences, and the LLR test in textbook form.

**Open / not covered here:**
- **The actual proofs.** Reproducing BLS Theorem 14, the Gerbicz-Li soundness bound, or the IBDWT error bound is out of scope — this doc cites them. Anyone *modifying* the algorithms (not just the code) must work from the references, not this summary.
- **PRST's exact parameter choices and their justification.** Why `genProthBase` picks the witness it does, the precise Kronecker-symbol `P` search, the BLS bound variants actually used, and the `Gerbicz_params` `log2b`-forced-to-1 simplification (`exp.cpp:390`) are engineering decisions whose *optimality* (vs. *correctness*) isn't analyzed here.
- **The smooth-exponent number theory.** *That* the smooth path computes `a^(b^n)` is documented (`exponentiation-algorithms.md` §1); *why* that specific quantity is the right thing for each test (how `k`, the tail correction, and `c` fold into the exponent) follows from the criteria above but isn't worked through candidate-by-candidate here.
- **Roots-of-unity attack algebra.** The specific forgeries `RootsTest` constructs and why the roots-of-unity check catches them is sketched in `test-harness.md` §4 and `proof-system.md` §7; the algebra is eprint 2023/195's.
