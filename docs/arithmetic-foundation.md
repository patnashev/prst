# Arithmetic foundation — deep dive

This is the layer everything else rests on. Two data types carry every number in PRST:

- **`Giant`** — an arbitrary-precision integer (the `k`, `b`, `c`, factor lists, residues, RES64s — all the *exact* integers). Backed by a swappable arithmetic strategy: a native giants library, GMP, or a per-GWnum-state allocator.
- **`GWNum`** — a number in GWnum's FFT/IBDWT representation, the form the *heavy modular multiplication* happens in. Created and operated on by a `GWArithmetic`, configured by a `GWState`.

The bridge between them — and the thing checkpoints persist — is **`SerializedGWNum`** (a `GWNum` flattened to portable words). And the layer that makes long FFT computations trustworthy is **`ReliableGWArithmetic`**, which detects round-off error and escalates. Every other doc in this folder has been quietly leaning on these: `Giant` is the type threaded through `inputnum-parsing.md`, the residues in `run-hierarchy.md`, and the records in `state-serialization.md`; `GWState` is the "math runtime" `main()` hands to every test.

This doc covers PRST's *own* arithmetic source (the `patnashev/arithmetic` submodule). The actual FFT bignum — **GWnum** (Woltman's `gwnum/`) and **GMP** — are prebuilt third-party libraries: PRST calls their C APIs (`gwsetup`, `gwmul3`, `mpz_*`) but owns no source, so they're understood at the interface here, not dissected. The wire (de)serialization, by contrast — `gwserialize`/`gwdeserialize`/`gwconvert` — is in-tree, defined in `arithmetic.cpp` (§5).

Source files:
- `framework/arithmetic/giant.{h,cpp}` (`Giant`, `GiantsArithmetic` + the `GW`/`GMP` backends, including the `GMPArithmetic` GMP binding)
- `framework/arithmetic/arithmetic.{h,cpp}` (`GWState`, `GWArithmetic`, `CarefulGWArithmetic`, `ReliableGWArithmetic`, `GWNum`, `SerializedGWNum`)
- `framework/arithmetic/integer.{h,cpp}` (small-integer number theory — `gcd`/`inv`/`is_prime`/`phi`/`jacobi`/`kronecker` — and the prime sieve `PrimeList`/`PrimeIterator`)

Prereqs / companions: `inputnum-parsing.md` (`InputNum::value()` returns a `Giant`; `InputNum::setup` drives `GWState::setup` and supplies `known_factors`), `state-serialization.md` (the `Giant` wire encoding and `SerializedGWNum` are what `Writer`/`Reader` move), `task-lifecycle.md` (the `ReliableGWArithmetic` reliable-mode contract it treats as a black box), `run-hierarchy.md` / `proof-system.md` (every test computes on `GWNum`s).

## 1. Data model

**`Giant` on the heap.** A `Giant` *is* a `giant_struct` (`giant.h:169-176`) wrapped in field-element machinery:

```cpp
struct giant_struct {
    int _capacity = 0;        // allocated uint32 words
    int _size = 0;            // word count; SIGN of _size = sign of the number
    uint32_t* _data = nullptr; // little-endian base-2^32 limbs; _data[0] is least significant
};
class Giant : public FieldElement<GiantsArithmetic, Giant>, protected giant_struct { … };
```

The field order is not incidental. `giant.cpp:53` defines `#define giant(x) ((::giant)&(x)._capacity)` — it reinterprets the address of `_capacity` as a `::giant` (gwnum's C bignum pointer) and hands it to the native ops (`itog`, `gianttogw`, `mulg`, …), which then read/write the struct in place. For that cast to be valid, `giant_struct`'s layout from `_capacity` onward must match whatever `::giant` points to — so the field order is load-bearing; reordering or inserting a field before them silently corrupts every native call. (The exact C `giant` definition lives in gwnum's `giants.h` — external/prebuilt — so the match is enforced by the cast, not visible here.) (`size()`/`capacity()` translate to 32-bit-word counts; under GMP the stored limbs are `GMP_NUMB_BITS` wide and `size()` converts + trims leading zeros — `giant.cpp:207-228`.)

**The arithmetic strategy.** `Giant` carries no logic; every operation delegates to its `GiantsArithmetic& arithmetic()`. There are four backends (`giant.h`):

| Backend | Role |
|---|---|
| `GiantsArithmetic` | the native giants library (default fallback) |
| `GMPArithmetic` (`giant.{h,cpp}`, `#ifdef GMP`) | GMP-backed (the `mpz_*` binding); **the default when `libgmp` loads** (`giant.cpp:23-37`) |
| `GWGiantsArithmetic` | native, but allocating out of a `gwhandle`'s memory (capacity sized to the FFT) |
| `GWGMPArithmetic` | GMP variant of the above |

`GiantsArithmetic::default_arithmetic()` is the process-wide singleton a bare `Giant x;` uses; `alloc_gwgiants(gwdata, capacity)` makes the per-`GWState` allocator (`giant.cpp:44-51`). So "which arithmetic" is chosen at construction and threaded through every op — a `Giant` and its arithmetic are inseparable.

**`GWNum` and friends.** A `GWNum` (`arithmetic.h:225`) wraps a single `gwnum _gwnum` (an FFT-domain buffer, `gwalloc`/`gwfree`). It's operated on by a `GWArithmetic` (the math) over a `GWState` (the config + handle). `SerializedGWNum` (`arithmetic.h:207`) is a `std::vector<uint32_t>` holding the `gwserialize` output of a `GWNum` — the only `GWNum` form that survives to disk. `GWNumWrapper` (`arithmetic.h:345`) is a non-owning view (doesn't free on destruct).

## 2. `GWState::setup` — the central method

`GWState` is the math runtime: it owns the `gwhandle handle`, the modulus `N`, the per-state `giants` allocator, the `fingerprint`, and the FFT description. The most illuminating method is the `k·b^n+c` setup (`arithmetic.cpp:58-84`), verbatim:

```cpp
void GWState::setup(uint64_t k, uint64_t b, uint64_t n, int64_t c)
{
    if (k >= (1ULL << 51) || b >= (1ULL << 32) || n >= (1ULL << 32) || std::abs(c) >= (1ULL << 30))
        throw ArithmeticException();                       // the fast-form domain limits
    Giant tmp;
    tmp.arithmetic().init((uint32_t*)&k, 2, tmp);
    tmp = tmp*power((Giant() = (uint32_t)b), (uint32_t)n) + c;   // N = k*b^n + c, as a Giant
    init();
    if (gwsetup(gwdata(), (double)k, (uint32_t)b, (uint32_t)n, (int32_t)c))   // configure the IBDWT FFT
        throw ArithmeticException();
    bit_length = (int)gwdata()->bit_length;
    if (gwdata()->GENERAL_MOD)
        bit_length /= 2;
    giants.reset(GiantsArithmetic::alloc_gwgiants(gwdata(), (bit_length >> 5) + 10));  // FFT-sized Giant allocator
    N.reset(new Giant(std::move(tmp)));
    if (gwdata()->GENERAL_MOD)
        bit_length = N->bitlen();
    fingerprint = *N%3417905339UL;                         // SAME modulus as InputNum::fingerprint
    if (!known_factors.empty() && known_factors > 1 && *N%known_factors != 0)
        throw ArithmeticException();
    if (!known_factors.empty() && known_factors > 1)
        *N /= known_factors;                               // divide out the known algebraic cofactor
    char buf[200];
    gwfft_description(gwdata(), buf);
    fft_description = buf;                                  // e.g. "FMA3 FFT length 512K"
    fft_length = (int)gwfftlen(gwdata());
}
```

Three things to read out of this:
1. **The fast-form domain is bounded** — `k < 2^51`, `b`/`n < 2^32`, `|c| < 2^30`. Outside it, callers use `setup(const Giant&)` (`gwsetup_general_mod`, a slower generic modulus) or `setup(int bitlen)` (`gwsetup_without_mod`, no modulus at all). `InputNum::setup` chooses among these (`inputnum-parsing.md` §4).
2. **`fingerprint = *N % 3417905339`** is computed here from the *actual* modulus GWnum built — and it's the **same constant** `InputNum::fingerprint()` uses. That's the cross-check: if the FFT was configured for the wrong number, the fingerprints disagree and the test aborts (the `ArithmeticException` you saw at the tail of `InputNum::setup`).
3. **`known_factors` is the algebraic-cofactor trick made concrete** — when `InputNum` recognized a `Phi`/`Quad`/`Hex` form, it set up GWnum modulo the *larger* fast-form multiple and passed the cofactor as `known_factors`; here it's divided out (`*N /= known_factors`), leaving the candidate as the working modulus while keeping the fast FFT (`inputnum-parsing.md` §6).

## 3. Field & method reference

**`Giant`** (`giant.h:176-550`):

| Member | Note |
|---|---|
| `data()` / `size()` / `capacity()` | raw limbs; word count (sign-stripped, GMP-translated); allocated words. |
| `empty()` | `_data == nullptr` (uninitialized — distinct from value 0). |
| `bitlen()` / `bit(i)` | bit length / bit test. |
| `to_string()` / `to_res64()` / `digits()` | decimal; low-64-bits hex (`"%08X%08X"`, `giant.cpp:185`); digit count. |
| operators `+ - * / % << >> == < …` | each constructs a result and delegates to `arithmetic()`. Mixed-type overloads for `uint32_t`/`int64_t`/`uint64_t` via the `GIANT_INIT_ADD_SUB` macro. |
| `power` / `gcd` / `inv` / `powermod` / `kronecker` / `rnd` / `substr` | number-theory ops, all on the backend. |
| `= GWNum` / `to_GWNum` | convert to/from the FFT form (`gianttogw`; negatives get `N` added first — mod-`N` representation, `giant.cpp:195-205`). |

**`GWState`** (`arithmetic.h:13-80`) — config knobs (set before `setup`) and live state (set by `setup`):

| Field | Kind | Meaning |
|---|---|---|
| `thread_count`, `spin_threads`, `instructions`, `next_fft_count`, `safety_margin`, `force_mod_type`, `large_pages` | config | FFT engine selection / sizing / threading. |
| `maxmulbyconst` | config | upper bound for `setmulbyconst` (the Fermat base `a`); set by `Run::create`. |
| `known_factors` | config | algebraic cofactor to divide out (§2). |
| `information_only` | config | report the FFT and abort (the `-fft info` path). |
| `handle` (`gwhandle`) | live | the GWnum library state. |
| `N` | live | the modulus `Giant` (null after `setup(bitlen)`). |
| `giants` | live | the FFT-sized `Giant` allocator. |
| `fingerprint`, `fft_description`, `fft_length`, `bit_length` | live | identity + FFT facts. |
| `mod_gwstate` | live | the secondary `GWState` that backs `mod()` (the known-factor reduction). Lazily created on the first `mod()` call (`force_mod_type=2`; set up over `*N`, or `square(*N)` under `GENERAL_MOD`) and used on *every* `mod()` call; the internal `known_factors > *N` test only adds one extra nested `mod_gwstate->mod()` reduction pass (`arithmetic.cpp:153-176`). `mod()` runs when `need_mod()` is true, i.e. `known_factors > 1`. |

Methods: `setup(k,b,n,c)` / `setup(Giant)` / `setup(bitlen)` (the three FFT configs), `copy(state)` (config only — **not** the live handle/N; you must `setup` after — this is what `-batch` uses per candidate), `clone(state)` (the deeper path that `gwclone`s the live handle, reached via the `GWState(GWState&)` copy-constructor; PRST's own code has no direct caller today), `done()` (`gwdone`+`gwinit`, the teardown), `need_mod()`/`mod()` (known-factor reduction), `ops()` (FFT-count → normalized op count).

**`GWArithmetic`** (`arithmetic.h:87-152`) — the FFT math over a `GWState`: `add`/`sub`/`mul`/`square`/`div`/`inv` plus fused `addmul`/`muladd`/`mulsub`/`mulmuladd`/… each taking GWnum `options` (`GWMUL_*` flags: which operands to preserve, whether to normalize). `setmulbyconst`/`setaddin`/`setpostaddin` tune the next multiply (asserting `|a| ≤ maxmulbyconst`). `popg()` mints a `Giant` from the state's allocator; `N()` is the modulus; `carefully()` returns the careful variant. Subclasses:
- **`CarefulGWArithmetic`** — `gwmul3_carefully` etc.: exact, slower, round-off-free. The `CarefulExp` tasks use it for short/exact legs.
- **`ReliableGWArithmetic`** — round-off-checked with escalation (§5).

**`GWNum`** (`arithmetic.h:225-343`) — `= Giant` (`gianttogw`), `= SerializedGWNum` (`gwdeserialize`), `to_string`. Division is `a * b.inv()` (GWnum has no native divide). `*num` yields the raw `gwnum`.

## 4. Lifecycle: how a test uses the layer

A test's path through this layer (drives the cadence in `task-lifecycle.md`):

1. **Configure.** `main()` fills a `GWState`'s config fields from the CLI; `InputNum::setup(gwstate)` picks `setup(k,b,n,c)` / `setup(Giant)` / `setup(bitlen)` and verifies the fingerprint. Now `gwstate.N`, `giants`, `fft_description` are live.
2. **Allocate.** The test makes a `GWArithmetic gw(gwstate)` and `GWNum`s from it (`gwalloc`).
3. **Compute.** The exponentiation/Lucas loop calls `gw.mul(...)` (or the reliable/careful variant) millions of times — this is the wall-clock. `setmulbyconst(a)` folds the Fermat base into the squaring.
4. **Checkpoint.** On the `DISK_WRITE_TIME` cadence, the current `GWNum` is captured into a `SerializedGWNum` (and exact values into `Giant`s) inside a `TaskState`, written by `File::write` (`state-serialization.md`). On resume, `= SerializedGWNum` deserializes it back.
5. **Finish + reduce.** The result `GWNum` is converted to a `Giant` (`= GWNum`), reduced mod `N` if `need_mod()`, and compared / `to_res64`'d for the result line.
6. **Teardown.** `gwstate.done()`; for `-batch`, a fresh `clone`/`copy`+`setup` per candidate.

`copy()` vs `clone()` is a real distinction: `copy()` brings only the config fields (not the live handle/N), so batch (`batch.cpp:310`) does `gwstate_cur.copy(gwstate)` then `input.setup(gwstate_cur)` to build a fresh FFT per candidate (BOINC's resume likewise relies on `setup` recomputing everything). `clone()` is the deeper `gwclone`-the-handle path, reached only via the `GWState` copy-constructor — not called directly in PRST's own code today.

## 5. SerializedGWNum and ReliableGWArithmetic

**`SerializedGWNum`** (`arithmetic.cpp:878-895`) is the checkpoint bridge. Assigning from a `GWNum` calls `gwserialize` twice — once with a null buffer to learn the length (it returns `GWERROR_MALLOC` with `-len`), then again to fill a resized `vector<uint32_t>`:

```cpp
SerializedGWNum& SerializedGWNum::operator = (const GWNum& a)
{
    int len;
    if (gwserialize(a.arithmetic().gwdata(), *a, NULL, 0, &len) != GWERROR_MALLOC) throw ArithmeticException();
    _data.resize(-len);
    if (gwserialize(a.arithmetic().gwdata(), *a, _data.data(), (int)_data.size(), &len) != 0 || _data.size() != len) throw ArithmeticException();
    return *this;
}
```

`to_GWNum` is `gwdeserialize` (or `0` if empty). Both calls take the `gwhandle`, so a value is serialized/deserialized *through* a configured FFT. Unlike `gwsetup`/`gwmul3`, the serializer is **in-tree**: `gwserialize`/`gwdeserialize`/`gwconvert` are all defined right here in `arithmetic.cpp` (the ~470-line block at `arithmetic.cpp:898-1405`). The wire format is a 6-word header (`GWSERIALIZE_HEADER_SIZE`, `arithmetic.cpp:1033`) carrying flags (`GWSERIALIZE_FLAG_IRRATIONAL`/`GENERALMOD`/`MMGW_MOD`, `:1034-1036`) plus the stored `b`/`FFTLEN`/`k` and `NUM_B_PER_SMALL_WORD` (written at `arithmetic.cpp:1128-1132`), followed by the FFT words. Because the format records the source FFT's length and base, `gwdeserialize` handles the cross-FFT-length question itself — see §8.

**`ReliableGWArithmetic`** (`arithmetic.cpp:528-592`) is the answer to "FFT multiplication is approximate — how do we trust a billion-squaring result?" Its `mul` escalates in three tiers, tracking an op counter:

```cpp
void ReliableGWArithmetic::mul(GWNum& a, GWNum& b, GWNum& res, int options)
{
    bool suspect = _suspect_ops.count(_op) != 0;
    if (!suspect) {
        gwerror_checking(gwdata(), true);
        gwmul3(gwdata(), *a, *b, *res, options);             // 1. normal, error-checked
        gwerror_checking(gwdata(), false);
        if (gw_get_maxerr(gwdata()) > _max_roundoff) {        // round-off too high (> 0.4)
            gw_clear_maxerr(gwdata());
            _suspect_ops.insert(_op);
            if (*res == *a || *res == *b) _restart_flag = true;  // in-place → can't redo here
            else suspect = true;
        }
    }
    if (suspect) {
        GWNum s1 = a; GWNum s2 = b;
        gwmul3(gwdata(), *a, *b, *res, options);              // 2. retry (maybe a HW glitch)
        if (gw_get_maxerr(gwdata()) > _max_roundoff) {
            gwmul3_carefully(gwdata(), *s1, …, *res, …);       // 3. careful (exact) retry
            if (gw_get_maxerr(gwdata()) > _max_roundoff) _failure_flag = true;  // FFT too small
        }
        else _suspect_ops.erase(_op);                         // it was a HW glitch, not a real suspect
    }
    _op++;
}
```

The contract (the surface `task-lifecycle.md` §8 treats as a black box):
- **`_op`** — monotonic operation index; the position in the computation.
- **`_suspect_ops`** — op-indices that once exceeded round-off; on a re-run they're done carefully from the start.
- **`restart_flag()`** — an in-place op overflowed and can't be redone in place ⇒ the caller **must** restart from the last checkpoint (this time the suspect op is in the set). Ignoring it silently keeps a bad result.
- **`failure_flag()`** — even careful arithmetic overflowed ⇒ the FFT is too small; bump `next_fft_count` and re-setup.
- **`reset()`** clears everything; **`restart(op)`** rewinds `_op` to a checkpoint position but keeps `_suspect_ops`. `_max_roundoff = 0.4`.

## 6. Pitfalls

- **`Giant`'s field order is an ABI contract with gwnum's C `giant`.** The `giant(x)` cast (`giant.cpp:53`) reinterprets `&_capacity` as a `::giant*`. Reordering `_capacity`/`_size`/`_data`, or adding a field before them, silently corrupts every native giants call. They live where they do *because* the C struct does.
- **`_size` is signed; `size()` is not.** The sign of `_size` is the number's sign; `size()` returns the magnitude word count. Reading `data()` without consulting the sign mishandles negatives — and `to_GWNum` adds `N` to a negative `Giant` first (mod-`N` form). The on-disk encoding (`state-serialization.md`) carries the sign in the length field for the same reason.
- **`SerializedGWNum` is tied to its modulus, but cross-FFT-length loading is defined in-tree.** It's the `gwserialize`/`gwdeserialize` wire form (both in `arithmetic.cpp`), always handled through a `gwhandle`, and meaningful only for the *number* it was serialized for. Loading into a *different-length* FFT is handled explicitly by `gwdeserialize` (§8): same length is the fast path, a **smaller** target transform throws (`"Can't deserialize to smaller transform."`), and a **larger** one re-radix-converts. So a checkpoint survives a reliable-mode FFT *bump* (larger), which is exactly the direction the bump goes.
- **`restart_flag` / `failure_flag` are not optional.** They're the only signal that an FFT result is untrustworthy; the `Task` loop polls them and restarts / enlarges the FFT. New code driving `ReliableGWArithmetic` directly that forgets to check them will commit corrupt residues.
- **`GWState::copy` is config-only.** It copies thread count, FFT hints, `known_factors`, etc., but **not** `handle`/`N`/`giants`. A copied state is unusable until `setup()`. (`clone()` is the one that brings the live FFT across, via `gwclone`.)
- **GMP vs native changes `size()`/`capacity()` semantics.** Under GMP the stored limbs are `GMP_NUMB_BITS` wide; `size()`/`capacity()` translate to 32-bit words on the fly (`giant.cpp:207-228`). Code poking `data()` must know which backend is active — the default is GMP when `libgmp` loads.
- **`maxmulbyconst` gates `setmulbyconst`.** `setmulbyconst(a)` asserts `|a| ≤ state.maxmulbyconst` (`arithmetic.h:136`). `Run::create` sets `options.maxmulbyconst` to the base; setting a base above it trips a `GWASSERT`.

## 7. Quick reference

| You want to… | Call |
|---|---|
| A scratch big integer | `Giant g;` (uses `default_arithmetic`, GMP if present) |
| Configure the FFT for `k·b^n+c` | `gwstate.setup(k, b, n, c)` (bounds: k<2⁵¹, b,n<2³², \|c\|<2³⁰) |
| …for an arbitrary modulus | `gwstate.setup(someGiant)` (general mod, slower) |
| …with no modulus | `gwstate.setup(bitlen)` |
| Per-candidate FFT in a batch | `dst.copy(src)` then `input.setup(dst)` |
| Multiply in the FFT domain | `gw.mul(a, b, res, options)` (or `gw.square`) |
| Exact / round-off-free op | `gw.carefully().mul(...)` |
| Round-off-checked op + retry | `ReliableGWArithmetic`; **check `restart_flag()`/`failure_flag()`** |
| `GWNum` → exact integer | `Giant g = gwnum;` (adds `N` if negative) |
| Persist a `GWNum` | assign to a `SerializedGWNum` (loadable into a same-or-larger FFT; §8) |
| Low 64 bits as hex (RES64) | `g.to_res64()` |
| Identity / FFT facts | `gwstate.fingerprint`, `.fft_description`, `.fft_length` |

## 8. Open questions / non-coverage

- **GWnum itself.** `gwsetup`/`gwmul3`/`gwmul3_carefully`/`gwfft_description` and the FFT-size/instruction-set selection are Woltman's prebuilt library (`gwnum/`); PRST uses the API and the `gwhandle` fields (`GENERAL_MOD`, `bit_length`, `mulbyconst`) but doesn't own the source. Treated as the interface contract, not dissected. Likewise GMP via `giant.cpp`'s `GMPArithmetic` — the *binding* is in scope (and partly shown), the GMP internals aren't. (Note: the (de)serialization path — `gwserialize`/`gwdeserialize`/`gwconvert` — is *not* external; it's defined in this repo's `arithmetic.cpp`, see §5/§8.)
- **The `giant.cpp` GMP binding in full.** This doc covers `GMPArithmetic`'s *role* (default backend) and the `size()`/`capacity()` translation; the per-op GMP glue (`mpz_*` calls for `mul`/`gcd`/`powermod`/`kronecker`/…) is the bulk of `giant.cpp`'s mechanical binding, summarized not enumerated.
- **The wider arithmetic library.** Elliptic-curve (`edwards`/`montgomery`/`group`), polynomial (`poly`), and Lucas-sequence (`lucas`) arithmetic build *on* `Giant`/`GWNum` but are their own subsystem → `curves-and-polynomials.md`.
- **`SerializedGWNum` cross-FFT portability — resolved.** Can a value serialized under FFT length L be `gwdeserialize`d into a length-L′ FFT for the same modulus? `gwdeserialize` (`arithmetic.cpp:1141-1405`) answers it in-tree by comparing the header's stored `FFTLEN`/`b`/flags against the live `gwdata`: an **identical**-length match is the fast path that copies the FFT words directly (`arithmetic.cpp:1156-1181`); deserializing into a **smaller** transform throws `"Can't deserialize to smaller transform."` (`arithmetic.cpp:1183-1186`); a **larger** transform re-radix-converts word-by-word into the wider layout (`arithmetic.cpp:1188-1395`), with extra carefully-multiplied correction factors for the `GENERAL_MOD`/`MMGW_MOD` mode changes. So a checkpoint resumes cleanly across a reliable-mode FFT *bump* (which only grows the transform); the `task.cpp:140-146` continue-from-`_state` path relies on exactly the larger-transform branch.
- **Why the round-off threshold is 0.4, and the IBDWT error theory.** `_max_roundoff = 0.4` is a tuning constant; the numerical-analysis justification (why that bounds a correct convolution) is GWnum/multiplication-theory territory → `math-and-theorems.md` and the `mult_*.pdf`s.
- **`ThreadSafeGWArithmetic` / `PolyMult`.** Friend classes named in `arithmetic.h` (`GWNum`'s friends) but not part of the single-`Task` main path; relevant only if multi-`Task` parallelism is added (cross-ref the parked thread-safety question in lay-of-the-land).

## 9. `integer.cpp` — number theory & the prime sieve

`integer.{h,cpp}` is the small-integer companion to the big-`Giant` layer: 32/64-bit number theory plus a process-wide prime generator. Nothing here is GMP glue (that lives in `giant.cpp`, §1) — it's plain `uint32_t`/`uint64_t` math that the candidate-parsing and factor-walk code leans on.

**Free functions** (`integer.h:9-25`, `integer.cpp:8-162`). Header-declared `uint32_t` helpers with `int` overloads: `gcd`/`inv` (Euclid and modular inverse), `is_prime` and `phi` (both trial-divide via the iterator below — `is_prime` up to `√a`, `phi` factoring `N` as it goes), and the Legendre/Jacobi/Kronecker symbols `jacobi`/`kronecker` (`integer.cpp:66-162`, the classic bit-twiddling implementations that fall back to returning a `gcd` when the arguments aren't coprime). These are interface-level utilities — small, exact, and stateless.

**`PrimeList`** (`integer.h:29-49`) is a sieved table of primes below a bound. The one that matters is the process-wide 16-bit cache: `PrimeList::primes_16bit()` (`integer.h:43`) lazily builds and returns a singleton `PrimeList(65536)` — every prime below 2¹⁶ (6542 of them, the same constant that opens the `_31bit_prime_count` table). The constructor (`integer.cpp:166-186`) is an odds-only bit-sieve; `sieve_range(start, end, list)` (`integer.cpp:226-230`) sieves an arbitrary `int` window on top of it.

**`PrimeIterator`** (`integer.h:51-99`) is the forward primes-as-a-stream view. `PrimeIterator::get()` (`integer.h:69`) starts at 2 over the 16-bit cache; `+= offset`, `++`, `find(prime)`, and `max()` walk it. Past the cached 16-bit range it transparently extends to 31-bit primes by sieving 64K-wide segments on demand (`init_range` → `PrimeList::sieve_range`, `integer.cpp:240-251`), and `sieve_range(start, end, list)` (`integer.cpp:314-316`) exposes a segmented 64-bit sieve over the shared `sieve_range_t` template (`integer.cpp:188-224`). To make seeking cheap, segment boundaries are precomputed from a large static delta-table: `_31bit_prime_count[]` (`integer.cpp:319-1343`) is folded by `init_31bit_prime_pos()` (`integer.cpp:1345-1357`) into a cumulative prime-position index `_31bit_prime_pos`. That static table is the bulk of the file — it is why `integer.cpp` is ~124 KB (1358 lines) despite the actual logic being a few hundred lines; treat the table itself as opaque generated data.

**Who consumes them.** `InputNum` parsing/factoring uses `PrimeIterator` and `kronecker` (`inputnum.cpp:311`, `:488`, `:1572`, `:1589`, `:1691`, `:1941`). The Pocklington and Morrison certificate walks pick witness primes / Lucas parameters through `PrimeIterator::get()` and `kronecker` (`pocklington.cpp:148`, `:579-582`; `morrison.cpp:167`, `:190`, `:393`, `:466`, `:830`). Others (`order.cpp`, `fermat.cpp`, `abc_parser.cpp`, the EC `group.cpp`) draw on the same free functions and iterator. So although it sits below the FFT machinery, this module is a real dependency of the test setup — not optional plumbing.
