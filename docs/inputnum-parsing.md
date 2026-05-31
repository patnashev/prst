# InputNum parsing — deep dive

`InputNum` is the parsed candidate. `main()` hands it a string like `"30006!4-1"`, `"5*2^100+1"`, `"Phi(3,10^20)"`, or `"(2^257-1)/535006138814359"`, and `InputNum::parse` turns it into a structured form: a *type* (`KBNC` / `FACTORIAL` / `PRIMORIAL` / `GENERIC`), the components `k`, `b`, `n`, `c`, an optional divisor `F`, a partial factorization of `N∓1`, and a pretty-printed display string. **`Run::create` dispatches almost entirely off this structure** (`input.type()`, `input.c()`, `input.n()`, `input.factors().size()`, `input.is_half_factored()`), so "the wrong test was selected" bug reports almost always trace back here.

Two things make this module trickier than a parser has any right to be:

1. **Classification is not purely syntactic.** `parse` produces a *tentative* type — it already collapses the trivial `k*b^n+c` shape (`k=1, n=1, /1`) to `GENERIC` at `inputnum.cpp:959-963` — then `process()` re-derives the rest: it can absorb a `k` divisible by `b` into the exponent, and detect Fermat / cyclotomic / quadratic / hexic structure on a `KBNC` candidate (`process` never promotes a plain/`GENERIC` number; the algebraic-detection block at `inputnum.cpp:1059` is gated on `_type==KBNC`). The type a user sees reported is the *post-`process`* type, not the one the tokens suggested.
2. **Factorials and primorials are stored as negative `Giant`s.** A factor list entry with `factor.first < 0` is a *pseudo-factor*: a packed encoding of `n!` or `n#` that `expand_factors()` later expands into real primes. Code that walks `_factors` and forgets to check the sign will multiply garbage.

Source files:
- `framework/inputnum.h` (the class, type constants, accessors)
- `framework/inputnum.cpp` (the parser, `process`, `setup`, `mod`/fingerprint, `factorize`, `expand_factors`)

Prereqs / companions: `run-hierarchy.md` (the consumer — `Run::create` keys off `type()`/`c()`/`factors()`/`is_half_factored()`/`expand_factors()`), and `state-serialization.md` (future; `read`/`write` use the `File` `Reader`/`Writer`).

## 1. The data model

Type constants (`inputnum.h:13-22`):

```cpp
static const int ZERO = 0;
static const int GENERIC = 1;
static const int KBNC = 2;
static const int FACTORIAL = 3;
static const int PRIMORIAL = 4;

static const int ALGEBRAIC_SIMPLE = 0;
static const int ALGEBRAIC_CYCLOTOMIC = 1; //     X^3 -+ 1 = (X^2 +- X + 1)(X -+ 1)
static const int ALGEBRAIC_QUAD = 2;       //  1/4 X^4 + 1 = (1/2 X^2 + X + 1)(1/2 X^2 - X + 1)
static const int ALGEBRAIC_HEX = 3;        // 1/27 x^6 + 1 = (1/3 X^2 + X + 1)(1/3 X^2 - X + 1)(1/3 X^2 + 1)
```

The private state (`inputnum.h:94-114`), with what each field holds:

| Field | Meaning |
|---|---|
| `_type` | `ZERO`/`GENERIC`/`KBNC`/`FACTORIAL`/`PRIMORIAL`. The *raw* type; `type()` overrides it to `GENERIC` when a divisor `_gf ≠ 1` is present. |
| `_gk`, `_gb` | `k` and `b` as `Giant`s (arbitrary precision). For `FACTORIAL`/`PRIMORIAL`, `_gb` is the *computed* `n!`/`n#` value. |
| `_n` | the exponent (or the factorial/primorial argument). |
| `_c` | the additive constant (`int64_t`; `\|c\|` capped at 63 bits during parse). |
| `_gf` | the divisor `F` from a `(<expr>)/F` cofactor form; `1` when absent. |
| `_factors`, `_cofactor` | the partial factorization of the quantity the deterministic test needs (`N-1` for `c=+1`, `N+1` for `c=-1`), plus the unfactored remainder. **Entries may be negative pseudo-factors** (see §6). |
| `_b_factors`, `_b_cofactor` | the factorization of `b` itself (KBNC only). `_factors` is seeded from these times `_n`. |
| `_custom_k/_b/_d/_f` | display overrides — the literal text the user typed, so the echo reads `Phi(3,10^20)` rather than the expanded integer. |
| `_gfn` | if `N` is a Fermat number `b^(2^gfn)+1`, the exponent; else `0`. |
| `_multifactorial` | the multifactorial step `m` (`n!m`); `1` for an ordinary factorial. |
| `_algebraic_type`, `_algebraic_k` | the detected algebraic form and its small parameter `k`; drive the `setup()` known-factors trick (§6). |

`N` itself is never stored as a field; it's recomputed on demand by `value()` (`inputnum.h:70`):

```cpp
arithmetic::Giant value() { return _type == GENERIC ? _gb : _type != KBNC ? (_gk*_gb + _c)/_gf : (_gk*power(_gb, _n) + _c)/_gf; }
```

So `GENERIC` *is* `_gb`; `FACTORIAL`/`PRIMORIAL` are `(k·b + c)/F` (with `_gb` the precomputed factorial); `KBNC` is `(k·b^n + c)/F`.

## 2. `parse` — the central function

`InputNum::parse` (`inputnum.cpp:504-993`) is ~490 lines. It runs in three movements; the line ranges below are the map.

**(a) Tokenize** (`:523-543`). A `Token` (`:385-457`) is either an operator (`+ - * / ^ ! #`) or an `Expr` (`:271-383`) — a number, an `name(args)` function call, or a parenthesized `()` sub-expression. The whole string is lexed into a `tokens` vector; leftover non-space input after the optional closing quote is `"excess symbols"`.

**(b) Pick a top-level shape.** In priority order:

- **Cofactor-divisor `(<expr>)/F` → GENERIC** (`:546-639`). Only taken when the first token is a `()` expression *and* every following operator is `/` with an `Expr` operand (the look-ahead loop at `:549-562`). It recursively `parse`s the inner expression, multiplies all the divisors into `gf`, checks `recursive.value() % gf == 0` (else `"not divisible"`), copies the recursive result's fields, and **leaves `_gf` set** so `type()` will report `GENERIC`. This is why `(2^257-1)/53500…` is tested with a plain probabilistic Fermat test, not a Pocklington test.

  ```cpp
  if (it_expr->expr && it_expr->expr->str_func == "()")
  {
      bool bf = false;
      auto itF = it_oper;
      while (itF != tokens.end() && itF->oper == '/')   // every following op must be '/'
      {
          itF++;
          if (itF == tokens.end() || !itF->expr) { bf = false; break; }
          bf = true;
          itF++;
      }
      if (itF != tokens.end()) bf = false;              // ...and nothing else after
      if (bf) { /* recursive parse, accumulate gf, copy fields, keep _gf */ }
  }
  ```

- **`Phi(3,x)` / `Phi(6,x)`** (`:644-689`) and **`Quad(x)` / `Hex(x)`** (`:690-756`) — algebraic constructors. Each recursively parses its argument as `…+1`, derives `k`, sets `c = 1`, and records `_algebraic_type`/`_algebraic_k` when the small `k` is in range (`Phi`: `k<1000`; `Quad`: `k<100`; `Hex`: `k<32`). `Phi(6,x)` is `Phi(3,-x)`; a leading `-` in the argument flips 3↔6 (`:651-655`).

- **The `k*b^n+c` walk** (`:761-966`) — the common path. It is a hand-rolled left-to-right scan, not a precedence parser:
  1. A `*` chain accumulates the `k` factors (`:764-793`), merging each operand's factorization with power `+1`.
  2. The next `Expr` is `b`; an optional `^n` sets the exponent (`:802-833`); `b`'s factorization is merged into `_factors` with power `n`.
  3. `!` → `FACTORIAL` (with a following `!`-run or a second arg giving the multifactorial step) (`:835-866`); `#` → `PRIMORIAL` (`:867-885`). Both compute `_gb` via the `factorial()`/`primorial()` helpers and push a negative pseudo-factor.
  4. A `/d` run divides out cofactors, merging with power `−1` and checking divisibility (`:899-931`).
  5. A trailing `+c` / `-c` sets `c` (`:943-958`); `\|value\|` over 63 bits is `"C too big"`.
  6. **The GENERIC fallback** (`:959-963`): if there's no constant and the shape is the trivial `k=1, n=1, /1`, the type degrades to `GENERIC` (it's just a bare number or expression):
     ```cpp
     else if (it_oper == tokens.end() && type == KBNC && gk == 1 && n == 1 && gd == 1)
     {
         type = GENERIC;
         n = 0;
     }
     else if (c_required)
         return InputNum::ParseResult(false, …, "C expected");
     ```

**(c) Commit + `process`** (`:972-992`). The locals are moved into the members, `_factors` is sorted (positive-power factors first, then by value, `:989`), and `process(gd, true)` finishes classification. `parse` returns a `ParseResult` (`inputnum.h:24-35`) carrying `success`, an error `pos`, and a `message` — the `pos` is what lets `main()` point a caret at the offending character.

Note the `c_required` parameter: top-level parses require a trailing `±c` (default `true`). The recursive parses for the `p(index)` function (`:298`) and `(…)` sub-expressions (`:321`) pass `c_required = false`, so a bare expression is legal there. `Phi`/`Quad`/`Hex` don't use that flag — they *append* `"+1"` to the argument before recursing (`:657, :703`), guaranteeing a constant; and the cofactor-divisor inner parse (`:566`) uses the default, relying on the user's inner expression to carry its own `±c`.

## 3. Field & method reference

The accessors are where the non-obvious overrides live. Every public method (`inputnum.h:56-87`):

| Method | Returns / note |
|---|---|
| `empty()` | `_gb == 0` — nothing parsed yet. |
| `type()` | `_gf == 1 ? _type : GENERIC` — **a divisor form always reports `GENERIC`**, regardless of the inner shape. |
| `k()` | `k` as `uint64_t`, or **`0` if it doesn't fit in ≤52 bits / 2 limbs**. Callers that test `k() == 1` must know a huge `k` reads as `0`. |
| `b()` | `b` as `uint32_t`, or `0` if `>1` limb. (`Run::create`'s Proth check tests `b() == 2`; a giant base reads as `0`.) |
| `n()` | `_n`, but `1` when `GENERIC`. |
| `c()` | `_c`, but `0` when `_gf ≠ 1`. |
| `f()` | `F` as `uint32_t` (0 if multi-limb). |
| `gk()`, `gb()`, `gf()` | the full-precision `Giant&`s. |
| `gfn()` | Fermat-number exponent, `0` if not a Fermat number. |
| `multifactorial()` | the `!m` step. |
| `algebraic_type()`, `algebraic_k()` | the detected algebraic form. |
| `value()` | recomputes `N` (see §1). Allocates — don't call it in a hot loop. |
| `mod(m)` | `N mod m` computed from the factor structure without materializing `N` (§5). |
| `fingerprint()` | `mod(3417905339)` — the 32-bit checksum (§5). |
| `bitlen()` | bit length of `N`, computed cheaply from `k`/`b`/`n`/`F` and only materialized when `≤ 64` (`:1390-1404`). |
| `is_half_factored()` | the BLS gate `Run::create` uses (§5). |
| `add_factor(g)` | move a factor from a cofactor into the factored list (used when an external sieve supplies a factor). |
| `b_factors()`, `b_cofactor()`, `factors()`, `cofactor()` | the factorization views. `factors()`/`cofactor()` are what `Run::create` and the deterministic tests read. |
| `factorize_minus1(depth)`, `factorize_small()` | trial-division helpers returning small divisor lists. `factorize_small` feeds trial division (`prst.cpp:295,312`) and `-batch` (`batch.cpp:261,278`); `factorize_minus1` feeds the proof security seed (`proof.cpp:104`). |
| `expand_factors()` | expand factorial/primorial pseudo-factors into primes (§6). |
| `input_text()`, `display_text()` | the echoed forms; `display_text()` is elided to ~30 chars. |

## 4. Lifecycle: parse → process → setup

The pipeline a candidate flows through:

1. **`parse(s)`** → tentative type, components, partial factorization, then calls…
2. **`process(gd, factored)`** (`inputnum.cpp:1002-1217`) — the classifier. It:
   - (if not already factored) factorizes `_gk` and, for `KBNC`, `_gb`, seeding `_factors` from `_b_factors × _n` (`:1013-1056`); folds the divisor `gd` in with negative powers.
   - sorts `_factors` and rolls `_b_cofactor^n` into `_cofactor` (`:1037-1041`).
   - **detects algebraic structure** for `KBNC, c=1, gd=1, ALGEBRAIC_SIMPLE` (`:1059-1108`): a Fermat-number exponent `_gfn` when `k=1` and `n` is a power of two (`:1061-1062`); then, when `\|bitlen(k) − log2(b)·n\|` is small, it computes `b^n − k` and tags `ALGEBRAIC_CYCLOTOMIC` (`±1`), or `b^n/2 − k` → `ALGEBRAIC_QUAD`, or `b^n/3 − k` → `ALGEBRAIC_HEX`.
   - normalizes the representation: reduces a common GCD of `b`'s exponents into `n` (`:1115-1132`), absorbs a `k` divisible by `b` into `n` (`:1134-1139`), and works the divisor `gd` back through `k`/`b` (`:1142-1214`).
   - rebuilds `_input_text` and `_display_text` via `build_text` (`:1219-1306`).
3. **`setup(GWState&)`** (`inputnum.cpp:1308-1388`) — hands the number to GWnum. For an algebraic form it sets `state.known_factors` to the algebraic cofactor and configures GWnum modulo a *fast-form multiple* of the candidate (e.g. cyclotomic uses `k³·b^(3n) ± 1`, `:1320-1325`). GWnum then divides the known factor out (`arithmetic.cpp:76-79`: `*N /= known_factors`), leaving the candidate as the working modulus. The payoff is **FFT speed**, not size: GWnum's irrational-base transform is fastest on the special `k·b^n±c` shape, and the fast-form multiple is actually *larger* than the candidate (§6). After setup it verifies `state.fingerprint == fingerprint()` (or `*state.N % 3417905339 == fingerprint()`), throwing `ArithmeticException` on mismatch.

`read`/`write` (`inputnum.cpp:21-47`) are the checkpoint path: `read` always produces a `GENERIC` number from the stored `Giant`; `write` stores the numerator `k·b^n + c` directly (note: *not* divided by `_gf` — for the common `_gf == 1` case that equals `value()`).

There's also a **recursion** dimension: `parse` constructs throw-away `InputNum`s for sub-expressions (`(…)`, `Phi`/`Quad`/`Hex` arguments, the `p(index)` prime function at `:289-315`). Each recursive `parse`/`process` runs the full pipeline on the fragment, then the parent splices the fields in.

## 5. `mod` / fingerprint / `is_half_factored`

**`mod(modulus)`** (`inputnum.cpp:1433-1501`) computes `N mod m` *without* building `N`. For `KBNC` it multiplies the residues of `_cofactor` and each factor `p^e` (binary exponentiation of the per-factor residue, with reductions to keep products in 64 bits), then applies `+c` and divides by `F` via a modular inverse (`f_inv`). If `F` isn't invertible mod `m` it falls back to `value() % m`.

**`fingerprint() = mod(3417905339)`** (`inputnum.h:72`) is a checksum against a fixed 32-bit modulus (`3417905339`, ≈ `0.8·2³²`). It's used purely as a *sanity check*: `setup` recomputes the residue of the modulus GWnum actually built and throws if it disagrees with `fingerprint()`. It also namespaces checkpoint filenames (`File::unique_fingerprint(input.fingerprint(), …)` throughout the test classes). It is **not** cryptographic and not collision-free — see Pitfalls.

**`is_half_factored()`** (`inputnum.cpp:1503-1519`) is the BLS "more than half of `N∓1` is factored" gate that `Run::create` calls before committing to a deterministic Pocklington/Morrison test:

```cpp
bool InputNum::is_half_factored()
{
    if (std::abs(c()) != 1)
        return false;
    if (_cofactor.empty() || _cofactor.bitlen()*2 + 10 < bitlen())
        return true;
    if (_cofactor.bitlen()*2 > bitlen() + 10)
        return false;
    Giant tmp;
    tmp = 1;
    for (auto& factor : _factors)
        tmp *= power(factor.first, factor.second);
    if (_c == 1)
        return tmp > _cofactor;
    tmp -= 1;
    return tmp > _cofactor;
}
```

The cheap bit-length bounds short-circuit the common cases; only the ambiguous middle band materializes the factored product. (Note this multiplies `_factors` directly — it relies on the pseudo-factors having already been expanded by `Run::create`'s `expand_factors()` call, since this runs *after* it on the deterministic path.)

## 6. Algebraic factoring & the negative-`Giant` factor encoding

**Algebraic forms.** GWnum is fastest on the special shape `k·b^n±c` (an irrational-base discrete weighted transform). A candidate such as `Phi(3,b^n) = (b^(3n) − 1)/(b^n − 1)` isn't of that shape, but it *divides* one that is. PRST recognizes the form (explicitly via `Phi`/`Quad`/`Hex` syntax in `parse`, or implicitly via the `b^n − k ≈ ±1` test in `process`) and, in `setup`, configures GWnum modulo the fast-form multiple and supplies the known algebraic factor via `state.known_factors`, which GWnum divides out — leaving the candidate as the working modulus while keeping the fast FFT. The multiple is **larger** than the candidate; the win is FFT efficiency, not a smaller number. The four cases and the fast-form multiple each sets up (`inputnum.cpp:1318-1349`):

| `_algebraic_type` | Recognized form | GWnum fast-form multiple |
|---|---|---|
| `ALGEBRAIC_CYCLOTOMIC` | `Phi(3/6, …)`, or `b^n − k = ±1` | `k³·b^(3n) ± 1` |
| `ALGEBRAIC_QUAD` | `Quad(…)`, or `b^n/2 − k = ±1` | `k⁴·b^(4n)/4 + 1` (parity-adjusted) |
| `ALGEBRAIC_HEX` | `Hex(…)`, or `b^n/3 − k = ±1` | `k⁶·b^(6n)/27 + 1` (parity-adjusted) |

If GWnum falls back to a general modulus (`GENERAL_MOD`/`GENERAL_MMGW_MOD`) — i.e. the fast form bought nothing — `setup` abandons the trick, resets to `ALGEBRAIC_SIMPLE`, and re-runs on the full number (`:1352-1359`).

**Factorials and primorials as negative `Giant`s.** This is the encoding most likely to bite. When `parse` hits `!` or `#`, it builds a *negated* small `Giant` and pushes it as a factor:

```cpp
// factorial n!m  (inputnum.cpp:861-865)
Giant factor;
factor = ((uint64_t)multifactorial << 32) + n;  // limb0 = n, limb1 = m
factor.arithmetic().neg(factor, factor);         // store NEGATED
gb = factorial(factor);
::add_factor(factors, factor, 1);

// primorial n#  (inputnum.cpp:875-879)
Giant factor;
factor.arithmetic().neg(it_expr->expr->value, factor);  // negated n
gb = primorial(factor);
n = factor.data()[0];
::add_factor(factors, factor, 1);
```

So a factor-list entry with `factor.first < 0` is a *pseudo-factor*: a one-limb negated value is `n#` (primorial), a two-limb negated value is `n!m` (factorial, `data()[0]=n`, `data()[1]=m`). The `factorial()` / `primorial()` free functions (`:459-502`) decode and compute the actual product. **`expand_factors()`** (`:1545-1606`) is what turns these into real primes — it sieves every prime up to `n` (or every multifactorial term), accumulates exponents into a `map`, and rebuilds `_factors`. `Run::create` calls it right before the half-factored test, which is why the deterministic-test path sees expanded primes while the raw parse result still carries the compact pseudo-factor.

**The `factorize` sieve** (`inputnum.cpp:96-269`) is the workhorse behind both `process` and the `factorize_*` helpers. It pulls out the power of two, then runs a segmented bit-sieve of small primes up to roughly `2^(2s)` (default `s=10`), trial-dividing `N` and moving the unfactored remainder to `cofactor`. An optional `is_factor` predicate lets callers substitute a cheaper test — e.g. `factorize_minus1` passes `mod(p) == 1` so it can find factors of `N−1` without computing `N−1` (`:1529`).

## 7. Pitfalls

- **A `(<expr>)/F` form is always `GENERIC`.** `type()` returns `GENERIC` whenever `_gf ≠ 1` (`inputnum.h:57`), so `Run::create` routes every cofactor-divisor candidate to a probabilistic Fermat test — never Proth, Pocklington, or Morrison — no matter how nicely the inner expression factors. If a user expects a deterministic test on `(b^n−1)/F`, this is why they don't get one.
- **`k()` / `b()` / `f()` silently return `0` for multi-limb values.** They're convenience narrowings (`inputnum.h:58-62`). `Run::create`'s `input.b() == 2` and `log2(input.gk()) < input.n()` mix the narrow and wide accessors deliberately; new dispatch code must not assume `k()`/`b()` reflect the true magnitude.
- **Negative entries in `factors()` are pseudo-factors, not real factors.** Anything iterating `_factors` (a GCD product, a residue computation) must either run after `expand_factors()` or guard `factor.first < 0`. `mod`, `is_half_factored`, and `print_info` all assume expansion has happened on the deterministic path.
- **`type()` is post-`process`, not what the tokens said.** `parse` collapses the trivial `KBNC` shape to `GENERIC` (`inputnum.cpp:959-963`); `process` then can absorb `k` divisible by `b` into `n`, reduce a common exponent GCD, or tag an algebraic form (`process` itself never touches `_type`). Reason about the *reported* type, not the literal input.
- **`fingerprint()` is a 32-bit non-cryptographic checksum.** It guards against a GWnum mis-setup (wrong FFT, wrong modulus) and namespaces files; it is not collision-resistant. Don't treat a fingerprint match as proof two inputs are equal, and don't widen its role without revisiting every `File::unique_fingerprint` caller (a checkpoint-format concern — see `state-serialization.md`).
- **The algebraic `setup` trick mutates state and can self-disable.** On a GWnum general-modulus fallback, `setup` resets `_algebraic_type`/`_algebraic_k` to simple and recurses (`:1352-1359`). Code reading `algebraic_type()` after `setup` may see a different value than before it.

## 8. Quick reference

| Input syntax | Parsed `type()` | Notes |
|---|---|---|
| `123456789` | `GENERIC` | bare number; `value() == _gb` |
| `5*2^100+1` | `KBNC` | `k=5, b=2, n=100, c=1` (Proth-eligible) |
| `2^257-1` | `KBNC` | `k=1, b=2, n=257, c=-1` |
| `100!-1` / `100!2-1` | `FACTORIAL` | `_gb = 100!` (or `100!!`); negative pseudo-factor in `_factors` |
| `100#+1` | `PRIMORIAL` | `_gb = 100#`; negative pseudo-factor |
| `(2^257-1)/535006…` | `GENERIC` | divisor form; `_gf ≠ 1` ⇒ always `GENERIC` ⇒ probabilistic Fermat test |
| `Phi(3,10^20)` | `KBNC` | cyclotomic; `_algebraic_type = CYCLOTOMIC` |
| `Quad(10^20)` / `Hex(…)` | `KBNC` | `_algebraic_type = QUAD` / `HEX` |
| `p(1000000)` | (as recursed) | the 1,000,000th prime; capped at index 105097565 |

| You want to… | Call |
|---|---|
| Get the integer `N` | `value()` (allocates) |
| Cheaply test size | `bitlen()` |
| Compute `N mod m` without building `N` | `mod(m)` |
| Know which test will run | inspect `type()`, `c()`, `n()`, `factors().size()`, then `expand_factors()` + `is_half_factored()` — exactly what `Run::create` does |
| Expand `n!`/`n#` factors to primes | `expand_factors()` |
| Echo the user's input / a short form | `input_text()` / `display_text()` |

## 9. Open questions / non-coverage

- **The commented-out ECM/Fermat-probe block** in `factorize` (`inputnum.cpp:157-268`) sketches an Edwards-curve ECM stage and a base-3 Fermat probe on the cofactor. It's dead today; if a future version wants deeper automatic factoring of `b`'s cofactor, this is the seam. Not covered here beyond noting it exists.
- **`build_text` display logic** (`:1219-1306`) — the elision, the `_custom_*` overrides, the factorial/primorial/divisor formatting — is faithfully reproduced output formatting, not algorithm. Treated as a black box; read it directly if a display bug is reported.
- **The exact `factorize` sieve schedule** (segment sizing, the `s`/`j` growth at `:130-147`) is a performance-tuned trial-division bound. The contract is "small factors out, remainder in `cofactor`"; the schedule itself is a tuning parameter not audited here.
- **The GWnum `setup` math** (FFT selection, `known_factors`, `force_mod_type`) belongs to the GWnum/`GWState` layer; this doc covers only how `InputNum` *chooses* the modulus and verifies it via the fingerprint. See the GWnum notes referenced in lay-of-the-land's open questions.
- **`print_info`** (`:1635-1942`) is the `-info` diagnostic dump (small factors, algebraic factor display, N∓1 factor fractions, Kronecker symbols). It's reporting, not parsing; covered only in passing.
