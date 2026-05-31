# Config DSL — deep dive

Every PRST front-end — `main`, `batch_main`, `boinc_main`, `net_main` — declares its command line the same way: a fluent `Config` chain that *builds a tree*, then `parse_args(argc, argv)` walks `argv` against that tree, writing straight into `Options` / `GWState` / local variables. So to understand any entry point's flags you have to understand this one DSL. It also reads `.ini` files through the same tree.

There are **two parallel hierarchies** and the trick is keeping them apart:

1. **The runtime tree** — `ConfigObject` and its subclasses (`ConfigKey`, `ConfigKeyValue`, `ConfigGroup`, `ConfigExclusive`, `ConfigKeyList`). This is what actually exists at parse time and consumes `argv`.
2. **The builder layer** — the `*Setup` class templates (`ConfigGroupSetup`, `ConfigExclusiveSetup`, …), a CRTP fluent API whose only job is to *construct* the runtime tree as you chain calls. Once `parse_args` runs, the builder is gone; only the tree matters.

A reader who conflates the two gets lost: the verbs you type (`value_number`, `group`, `exclusive`, `end`) are builder methods; the objects they leave behind (`ConfigKeyValueNumber`, `ConfigGroup`, …) are what parse.

Source files:
- `framework/config.h` (both hierarchies — the runtime `ConfigObject`s and the `*Setup` builders)
- `framework/config.cpp` (`parse_args`/`parse_ini` per object type, `parse_number`)

Prereqs / companions: the entry-point docs that show real chains — `abc-batch.md` §3 (`batch_main`), `boinc-and-net.md` §2 (`boinc_main`), and `prst.cpp:82+` (`main`); `arithmetic-foundation.md` (`GWState`, a common parse target); `state-serialization.md` (`File`/`TextReader`, which `parse_ini` reads through).

## 1. The two hierarchies

**Runtime tree** (`config.h`), each node a `ConfigObject` with `parse_args(argc, argv, cur)` / `parse_ini(...)`:

| Node | Matches | Effect |
|---|---|---|
| `ConfigKey` (plain) | a bare flag token (`-batch`) | nothing — `ignore`. |
| `ConfigKeyCheck<T,V>` | a flag | `dst = value` (e.g. `-d` → `log_level = LEVEL_INFO`). |
| `ConfigKeyCode` | a flag | runs a `std::function<void()>`. |
| `ConfigKeyValue*` | `key`+`delim`+`value` | sets a target via `set_value`. Subtypes: `…Number<T,V>` (range-checked), `…String`, `…Enum<T,V>` (string→enum via `Enum`), `…Code` (`std::function<bool(const char*)>`). |
| `ConfigKeyList` | `key` then N values | fills N targets (one per registered value); `_list_delim` separates them (a char, or `' '` = successive `argv` tokens). |
| `ConfigGroup` | a name token, then greedily its sub-keys | a prefix namespace (`-fft`, `-check`, `-time`). |
| `ConfigExclusive` | the first of its `ConfigCase`s that matches | mutually-exclusive options; once one matches, only its `optional()` sub-keys are accepted. |
| `Config` | — | the root `ConfigGroup` + a `_default_code` for unmatched tokens. |

Every `ConfigObject` also carries an optional **`_parsed`** hook (a `ConfigKey`), run by `on_parsed()` *after* the object successfully parses — this is how `on_check`/`on_code` attach a side effect ("if `-fermat` appeared at all, set `ForceFermat = true`").

**Builder layer** (the `*Setup` CRTP templates). `ConfigGroupSetup<P, R>` carries the verbs (`value_number`, `check`, `group`, …); `P` is the concrete setup type each verb returns (so chaining stays typed), `R` is the parent returned by `end()`. `group()`/`list()`/`exclusive()` descend into child setups; `end()` climbs back. `Config` is *both* a `ConfigGroup` and a `ConfigGroupSetup<Config, nullptr_t>`, so the root is its own builder — that's why `Config cnfg; cnfg.value_number(...).group(...)…` works.

## 2. `parse_args` — the central methods

Two methods carry the dispatch. The **top loop** (`config.cpp:16-30`):

```cpp
void Config::parse_args(int argc, char *argv[])
{
    int cur = 1;
    while (cur < argc)
    {
        int last = cur;
        for (auto it = _objects.begin(); it != _objects.end(); it++)
            (*it)->parse_args(argc, argv, cur);     // each object advances cur iff it matches
        if (cur == last)                            // nothing matched argv[cur]
        {
            _default_code(argv[cur]);               // → the input number / batch file
            cur++;
        }
    }
}
```

Objects are tried **in registration order**; a match advances `cur`, an unmatched token falls to `_default_code`. The **matcher** that most flags use, `ConfigKeyValue::parse_args` (`config.cpp:224-256`), is where the `delim` convention lives:

```cpp
void ConfigKeyValue::parse_args(int argc, char *argv[], int& cur)
{
    if (cur >= argc) return;
    if (strncmp(argv[cur], _name.data(), _name.size()) != 0) return;  // argv must START WITH the key
    char* val_s = argv[cur] + _name.size();         // …the rest is the candidate value
    int inc = 1;
    if (_delim == ' ') {                             // value is the NEXT token:  "-t 4"
        if (*val_s != 0) return;                     //   key must be a whole token
        if (cur + 1 >= argc) return;
        val_s = argv[cur + 1]; inc = 2;
    }
    else if (_delim != 0) {                          // value follows a literal delim: "-fft+1"
        if (*val_s != _delim) return;
        val_s++;
    }
    else {                                           // glued:  "-t4"
        if (*val_s == 0) return;
    }
    if (!set_value(val_s)) return;                   // range/format check; failure = no match
    cur += inc;
    on_parsed();                                     // fire the on_check/on_code hook
}
```

Three `delim` modes, all in use:
- **`' '`** — value is the following `argv` token (`-t 4`, `-fermat a 3`).
- **`0`** — value is glued to the key in the same token (`-t4`).
- **any other char** — value follows that literal in the token (`-fft+1`, delim `'+'`).

`main`/`batch_main` register `-t` *both* ways (`value_number("-t", 0, …)` and `value_number("-t", ' ', …)`) so `-t4` and `-t 4` both work. A failed `set_value` (out of range, non-numeric, unknown enum) is treated as **no match**, so the token can fall through to another object — this is the main guard against prefix collisions (a `-t` value-key declines `-time`, because `set_value("ime")` fails the number parse). It's not foolproof: a longer key whose tail *is* a valid value could still be mis-claimed, so it works because the actual keys are chosen not to collide that way (see §6).

`ConfigKey::parse_args` (flags) is stricter: the rest after the name must be **empty** (`-d`, not `-debug`), then it runs `on_key()` (the side effect) and `on_parsed()`.

## 3. Builder-verb reference

All on `ConfigGroupSetup` (so available at the root and inside any `group`/`exclusive` case):

| Verb | Builds | Meaning |
|---|---|---|
| `ignore(key)` | `ConfigKey` | match a flag and do nothing (`-batch`/`-boinc` ignored by the sub-`main` that already dispatched on it). |
| `check(key, dst, value)` | `ConfigKeyCheck` | flag → `dst = value`. |
| `check_code(key, fn)` | `ConfigKeyCode` | flag → run `fn()`. |
| `value_number(key, delim, dst, min, max)` | `ConfigKeyValueNumber` | parse a number into `dst`, range-checked. `dst` may be `std::optional<T>`. |
| `value_string(key, delim, dst)` | `ConfigKeyValueString` | copy the raw string. |
| `value_enum(key, delim, dst, Enum)` | `ConfigKeyValueEnum` | map a string to an enum value via an `Enum<V>` (`.add("AVX","AVX")…`). |
| `value_code(key, delim, fn)` | `ConfigKeyValueCode` | hand the value to `fn(const char*) → bool` (e.g. `-order`, `-ini`). |
| `group(key)` | `ConfigGroup` → child setup | a prefix namespace; chain its sub-keys, then `end()`. |
| `list(key, delim, list_delim, fixed?)` | `ConfigKeyList` → child setup | a key followed by several values (e.g. `-proof name <pt> <prod>`); chain `value_*`, then `end()`. |
| `exclusive()` | `ConfigExclusive` → child setup | open a set of mutually-exclusive `ex_case()`s. |
| `ex_case()` / `optional()` | `ConfigCase` / its optional group | one alternative; `optional()` holds keys allowed only once the case matched. |
| `on_check(dst, value)` / `on_code(fn)` | sets `_parsed` on the **last-added object** | a side effect that fires when that object (often a whole `group`) parses. |
| `end()` | — | climb back to the parent setup. |
| `default_code(fn)` (root only) | sets `Config::_default_code` | the catch-all for unmatched tokens — the candidate expression or batch filename. |

`parse_number` (`config.cpp:190-208`) accepts an order suffix: `k`/`K`=10³, `M`=10⁶, `G`=10⁹, `T`=10¹², `P`=10¹⁵ — applied *before* the range check (so `-t 1k` parses 1000, then fails `max 256`).

## 4. How a chain parses

1. **Build.** The fluent chain runs at construction, leaving a `ConfigObject` tree under the root `Config`. Nothing parses yet.
2. **Walk.** `parse_args` loops `argv`; each pass tries every root object in order. A named **`ConfigGroup`** that matches its name token then runs its own fixpoint loop (`config.cpp:311-318`): re-scan its sub-objects until a full pass consumes nothing. So after `-fft`, the tokens `+1 safety 0.5 generic` can appear in any order and each is consumed by a sub-key; the first token that *no* sub-key matches ends the group.
3. **Exclusive.** A `ConfigExclusive` tries each `ConfigCase` until one matches; a case matches only if it consumes **all** its required objects (`ConfigCase::parse_args`, `config.cpp:471-499`). The winning case is latched in `_case`, and from then on only its `optional()` keys are accepted.
4. **Fall through.** Any token no object consumed in a full root pass goes to `_default_code`.
5. **Side effects.** Each successful parse calls `on_parsed()`, firing any `on_check`/`on_code` hook attached to that object.

**`.ini` files** (`parse_ini`, `config.cpp:32-87`) go through the same tree: `[section]` → group `-section`, `key=value` → a key under that group, bare `key` → a valueless flag, `;`/`#` → comment. Nested groups use `\` in section names. One normalization to know when writing an `.ini` (`config.cpp:57-59`): **top-level** keys are stored with a leading `-` prepended (so an `.ini` line `t=4` matches the `-t` flag), but keys **inside a `[section]`** are kept bare (matching the group's sub-key names like `safety`, `generic`). The tree is matched against the parsed `values` map; leftover entries print `Unknown option …`. (`-ini <file>` itself is a `value_code` that reads the file via `File`/`TextReader` and calls `parse_ini`.)

## 5. Common patterns

- **Dual `-t4` / `-t 4`.** Register the same key with `delim 0` and `delim ' '`. Order them glued-first; the glued one fails `set_value` on an empty value and falls through to the spaced one for `-t 4`.
- **Prefix group.** `group("-fft").value_number("+", 0, …).value_number("safety", ' ', …).check("generic", …).end()` — `-fft +1 safety 0.5 generic`. Note `"+"` with `delim 0` inside the group makes `+1` parse as the `+` key glued to `1`.
- **Group-level flag via `on_check`.** `group("-fermat").value_number("a", ' ', …).end().on_check(options.ForceFermat, true)` — `a` is optional, but the mere presence of `-fermat` sets `ForceFermat` (the hook fires on the group's own parse).
- **Exclusive switches.** `group("-check").exclusive().ex_case().check_code("near", […]).end().ex_case().check_code("always", […]).end()…end()` — `-check near|always|never`, first match wins.
- **Positional list.** `group("-proof").exclusive().ex_case().value_number("save", ' ', count, …).optional().list("name", ' ', ' ').value_string(pt).value_string(prod).end()…` — `-proof save N name <pt> <prod>`.
- **Catch-all.** `default_code([&](const char* p){ if (!input.parse(p)) printf("Unknown option %s.\n", p); })` — the candidate number; in `batch_main` it's the batch filename.

## 6. Pitfalls

- **Registration order is matching priority.** Objects are tried in the order added (`config.cpp:22`, `:315`); the first to consume a token wins that pass. Reordering a chain can change which key claims an ambiguous token. Value-keys are partly self-protecting (a failed `set_value` declines the token), but flags (`ConfigKey`, exact-match) and groups (exact name) are not order-sensitive only because they require the rest of the token to be empty.
- **`delim` is the most error-prone knob.** `0` = glued, `' '` = next token, anything else = literal separator. Picking `' '` when you meant `0` (or vice-versa) silently makes a flag never match (it waits for a separate token that never comes, or expects a glued value that isn't there).
- **`on_check`/`on_code` bind to the *last-added* object.** They must be placed immediately after the object (often a `group(...).end()`) they annotate; put elsewhere they attach to the wrong node. The hook fires on that object's *successful parse*, not unconditionally.
- **A partially-matched exclusive case still ran its side effects.** `ConfigCase::parse_args` commits `cur` only if it consumes *all* its objects, but the `set_value`/`on_key` calls of the ones it *did* match already mutated their targets — and (unlike `parse_ini`, which saves/restores copies, `config.cpp:533-538`) `parse_args` does **not** roll them back. Harmless for the usual single-key cases, latent for multi-key ones.
- **Order suffixes apply before range-checking.** `parse_number` multiplies by `k`/`M`/`G`/… first, so a suffixed value can overflow or exceed `max` and be *rejected* (declining the token) rather than clamped.
- **`set_value` failure is silent.** A bad value isn't an error; the key just doesn't match, and the token usually ends up in `_default_code` (or, for a bare number, mis-read as the candidate). A typo'd numeric flag can quietly become "the input."

## 7. Quick reference

`delim`: `0` glued (`-t4`) · `' '` next token (`-t 4`) · `c` literal (`-fft+1`). `set_value` fail ⇒ token declined (falls through).

| You want to… | Verb |
|---|---|
| A boolean/enum flag setting a field | `check(key, dst, value)` |
| A flag running code | `check_code(key, fn)` |
| A numeric option (range-checked, `k`/`M`/… suffixes) | `value_number(key, delim, dst, min, max)` |
| A string / enum option | `value_string` / `value_enum(…, Enum<V>().add(...))` |
| An option needing custom parsing | `value_code(key, delim, fn→bool)` |
| A prefix namespace of sub-keys | `group(key)…end()` |
| Several values after one key | `list(key, delim, list_delim, fixed?)…end()` |
| Mutually-exclusive alternatives | `exclusive().ex_case()…end()…end()` |
| "If this (group) appeared, also set X" | `…end().on_check(dst, value)` |
| The unmatched-token catch-all | `default_code(fn)` (root) |
| Read options from a file | `value_code("-ini", ' ', …parse_ini)` |
| Match-and-ignore a token | `ignore(key)` |

## 8. Open questions / non-coverage

- **The CRTP template plumbing.** The `P`/`R` parameters of `ConfigGroupSetup<P,R>` and the `ConfigKeyGroupSetup`/`ConfigExclusiveSetup`/`ConfigCaseSetup`/`ConfigOptionalSetup` derivations exist purely to keep the fluent chain's return types correct through `end()`. This doc explains the *behavior* (descend/climb); the exact template-instantiation chain is mechanical C++ not re-derived here.
- **`from_chars` portability.** `config.cpp` has two `from_chars` implementations gated on `CHARCONV` (`std::from_chars` vs. hand-rolled `strtoull`/`strtoll`/`strtod`); the numeric-parsing edge cases (`isnormal` rejection of floats, `ERANGE`) are noted but not exhaustively tabulated.
- **Exact INI section nesting.** `parse_ini`'s handling of `\`-separated nested section names and the `values[""]`/default-code interplay (`config.cpp:68-79`) is reproduced operationally; the full grammar of valid `.ini` (see `sample.ini`) isn't enumerated.
- **Per-flag semantics live with the caller, not here.** *What* `-fft +1` or `-check strong count N` means for the run is documented where the chain is declared (the entry-point docs) and where the target field is consumed (`arithmetic-foundation.md` for `GWState`, the test docs for `Options`); this doc covers only the parsing mechanism.
