# Checkpoints — PRST state types and files

What PRST persists between sessions: checkpoints (`.ckpt`), recovery points (`.rcpt`), progress params (`.param`), proof artifacts (`.proof.<i>`, `.cert`, `.pack`), and per-mode variants (`.div` children for `-divides`, `.ord<fingerprint>` for `-order`). Everything rides the framework's `Writer`/`Reader`/`File` framing — the binary header, atomic write, `.md5` sidecar, and the `TaskState` record shape are documented in the framework's `state-serialization.md`. This doc owns the PRST-specific layer: the `TYPE` registry, the state classes, the strong-check file handshake, and the LLR2 compatibility munging. The framework doc is deliberately consumer-agnostic — Prefactor, another consumer of the same framework, uses its own checkpoints, and they're completely different.

Source files:
- `src/exp.h` / `exp.cpp` (TYPEs 1, 2, 5, 8; the strong-check handshake)
- `src/lucasmul.h` (TYPEs 9, 10, 11)
- `src/proof.h` (TYPEs 3, 4, 6)
- `src/support.{h,cpp}` (`LLR2File`) and `src/net.cpp` (`LLR2NetFile`)
- `src/prst.cpp` (`File::FILE_APPID = 4`, the registry comment, the checkpoint filenames)

Prereqs / companions: the framework's `state-serialization.md` (the on-disk format and the `TaskState` base class), `task-lifecycle.md` (when checkpoints are written — the `on_state` cadence and `write_state`), `exponentiation-algorithms.md` (the tasks that own these states), `proof-system.md` (TYPEs 3/4/6 and the `.pack` container), `boinc-and-net.md` (the same buffers over HTTP).

## 1. The TYPE registry

The canonical registry is the comment at `src/prst.cpp:54-65`, kept next to `File::FILE_APPID = 4`. Mirrored here with the owning classes:

| TYPE | Meaning | Owner class | File | Body after `iteration` |
|---|---|---|---|---|
| −1 | test metadata | — (reserved file-level id, not a `TaskState`) | | — |
| 0 | number | — (reserved file-level id, not a `TaskState`) | | — |
| 1 | checkpoint | `BaseExp::StateValue` | `src/exp.h:42` | `Giant` value |
| 2 | strong check checkpoint | `StrongCheckMultipointExp::StrongCheckState` | `src/exp.h:308` | `int recovery` + `SerializedGWNum X` + `SerializedGWNum D` |
| 3 | proof product | `Proof::Product` | `src/proof.h:26` | `Giant X` (here `iteration` = tree depth) |
| 4 | certificate | `Proof::Certificate` | `src/proof.h:40` | `Giant X` + optional `Giant a_power` + `Giant a_base` |
| 5 | strong check placeholder | **no class** — set inline via `_state.reset(new TaskState(5))` | `src/exp.cpp:459`, `:572`, `:751`, `:778`; `src/lucasmul.cpp:220`, `:462`, `:473` | nothing — iteration only |
| 6 | proof checkpoint | `Proof::State` | `src/proof.h:59` | **no standard `read`/`write` override** (see notes) |
| 7 | — unused | a hole, not a free slot | | |
| 8 | serialized checkpoint | `BaseExp::StateSerialized` | `src/exp.h:29` | `SerializedGWNum` value |
| 9 | LucasV checkpoint | `LucasVMulFast::State` | `src/lucasmul.h:35` | `int index` + `Giant V` + `int parity` |
| 10 | LucasUV checkpoint | `LucasUVMul::State` | `src/lucasmul.h:110` | `Giant Vn` + `Giant Vn1` + `int parity` |
| 11 | LucasUV strong check checkpoint | `LucasUVMul::StrongCheckState` | `src/lucasmul.h:130` | `int recovery` + `SerializedGWNum Vn` + `Vn1` + `int Vparity` + `U` + `V` + `int parity` |

Notes:
- **TYPE 5 is an important state type.** It means the checkpoint is 0 iterations after the recovery point. Since it's empty, it does not have its own class — the record persists only the base-class iteration (§4 shows where it's installed).
- **TYPE 6 (`Proof::State`) is the one record without a `read`/`write` override** (`src/proof.h:56-77`). Through `File::read/write` it would persist only the base `iteration`; its `X`/`Y`/`exp`/`h` payload is managed by the proof code (`ProofSave`/`ProofBuild`). Don't assume the standard "iteration + fields" layout applies to it — see `proof-system.md`.
- `Proof::Certificate::read` is **forward/backward tolerant**: it reads `X`, then *optionally* `a_power`+`a_base` via `((reader.read(_a_power) && _a_power != 0 && reader.read(_a_base)) || true)` (`src/proof.h:48`) — a 1-field certificate for smooth numbers still loads. This is the one record that deliberately tolerates a shorter body; the rest fail closed on truncation.
- `version()` is `0` for every state today; `bool`s are written as `int` `1`/`0` (e.g. `parity`, `LucasVMulFast::State::write`, `src/lucasmul.h:45`).
- **Never reuse or renumber.** 7 is the only hole, and it is not a free slot; 5 is a live placeholder with no class of its own. New states append (≥ 12) and update the `prst.cpp` comment.

## 2. The state class tree

```
TaskState                                 (framework base; state-serialization.md)
├── (no class)                            TYPE=5   (strong check placeholder: iteration only)
├── BaseExp::State
│   ├── BaseExp::StateSerialized          TYPE=8   (mid-task: SerializedGWNum)
│   └── BaseExp::StateValue               TYPE=1   (final iteration: Giant)
├── StrongCheckMultipointExp::StrongCheckState   TYPE=2   (Gerbicz/Li check intermediates: _recovery, X, D)
├── LucasVMulFast::State                  TYPE=9   (single V plus index + parity)
├── LucasUVMul::State                     TYPE=10  (V_n, V_{n+1}, parity)
├── LucasUVMul::StrongCheckState          TYPE=11  (UV-form Gerbicz check intermediates)
├── Proof::Product                        TYPE=3   (proof-product checkpoint)
├── Proof::Certificate                    TYPE=4   (final certificate written by ProofBuild)
└── Proof::State                          TYPE=6   (proof checkpoint state)
```

The tasks that own these states (the `BaseExp` exponentiation family, the Lucas multipliers, the proof tasks) are toured in `exponentiation-algorithms.md` §1 and `proof-system.md`.

## 3. StateValue vs. StateSerialized, and `Point::value`

`BaseExp::State` can be a `StateValue` (TYPE 1, an exact `Giant`) or a `StateSerialized` (TYPE 8, a cheap FFT-domain `SerializedGWNum`), and which one a given checkpoint uses can be manipulated via the `MultipointExp::Point::value` flag. This is not related to file serialization or the on-disk format — it is PRST exponentiation logic, and it's documented where that logic lives: `exponentiation-algorithms.md` §1 (the `State` classes), §2 (`State::cast` picking the concrete type at each point), and its pitfalls (`commit_execute` chooses the type by iteration). The short version: point builders set `value = (pos == n)`, so only the final/selected points materialize a `Giant`; everything in between checkpoints the cheaper serialized form.

## 4. The `.ckpt`/`.rcpt` strong-check handshake

A strong-check task (Gerbicz / Gerbicz-Li) keeps **two** files, opened in `src/prst.cpp:403-404`:

- `prst_<fingerprint><suffix>.ckpt` — the working checkpoint: mid-block `StrongCheckState` (TYPE 2, or 11 for Lucas UV), or the TYPE 5 placeholder.
- `prst_<fingerprint><suffix>.rcpt` — the recovery point: the last *verified* `State` (TYPE 1/8), written only after a check passes.

`StrongCheckMultipointExp::init_state` (`src/exp.cpp:447-465`) installs the placeholder: when the current position *is* the recovery point — no unverified work in flight — the checkpoint slot gets `_state.reset(new TaskState(5)); _state->set(_state_recovery->iteration());` (`:459-460`). On resume, a TYPE 5 checkpoint tells the task "start exactly at the recovery point"; a TYPE 2 checkpoint carries the mid-block intermediates to continue from. A failed Gerbicz check rolls back by discarding the checkpoint and re-reading the recovery point — in-block progress is lost, verified progress never is.

`StrongCheckMultipointExp::write_state` (`src/exp.cpp:467-475`) writes the recovery file first (if the recovery state is unwritten), then delegates to `Task::write_state()` for the checkpoint. The same TYPE 5 reset appears mid-run whenever the recovery point catches up with the checkpoint (`src/exp.cpp:572`, `:751`, `:778`; Lucas analogues in `src/lucasmul.cpp`).

Sub-runs get child files: `file.add_child(name, File::unique_fingerprint(fp, salt))` scopes per-base/per-factor checkpoints (e.g. `-divides` creates `.div` children per base and power). The fingerprint math is the framework's (`state-serialization.md`); the salting conventions are PRST's.

## 5. LLR2 on-disk compatibility

`LLR2File` (`src/support.cpp:13-52`) provides bidirectional compatibility with LLR2's on-disk format, so a PRST worker can resume an LLR2 checkpoint and vice-versa. The munging is purely in the header/trailer bytes:

- **On read** (`:20-31`): if the buffer looks like an LLR2 file (`MAGIC_NUM`, **`_buffer[4] == 2`** — i.e. LLR2's `appid`), rewrite `_buffer[4] = 4` (PRST's appid) and `_buffer[6] = _type` (set the TYPE byte), and for a `StateValue` (TYPE 1) **decrement** the iteration word at offset 12.
- **On write** (`:33-52`): the inverse — `buffer[4] = 2`, `buffer[6] = 0`, **increment** the `StateValue` iteration at offset 12, then append an LLR2 trailer: a zero `uint32`, a 32-bit additive checksum (sum of `uint32`s from offset 8 to end), and 20 zero `uint32`s.

Byte 4 is the **appid**, not a format version — LLR2 stamps `2`, PRST stamps `4`. The off-by-one on the `StateValue` iteration reflects that LLR2 counts the iteration one differently from PRST.

`LLR2NetFile` in `src/net.cpp` is the same munging over an HTTP-backed `NetFile` (see `boinc-and-net.md` §5).

## 6. Pitfalls

- **The registry is a permanent contract.** Reusing or renumbering a TYPE silently aliases two on-disk formats; 7 is the only hole and not a free slot, 5 is live despite having no class. Append ≥ 12 and keep the `prst.cpp:54-65` comment in sync.
- **`.ckpt` and `.rcpt` are not interchangeable.** The checkpoint may hold unverified work; only the recovery point is check-verified. Deleting `.rcpt` and keeping `.ckpt` forfeits the rollback target (see the `exponentiation-algorithms.md` pitfalls for the in-memory analogue).
- **The LLR2 munging pokes fixed offset 12** — it assumes a fingerprinted file (body at offset 12). A fingerprint-0 file would put the iteration at offset 8; the LLR2 path never writes such files, but don't reuse the code for one.

## 7. Quick reference

| Artifact | File | TYPE(s) | Written by |
|---|---|---|---|
| Working checkpoint | `.ckpt` | 1, 8, 2, 9, 10, 11, or 5 (placeholder) | `Task::write_state` on the `on_state` cadence |
| Recovery point | `.rcpt` | 1, 8 | `StrongCheckMultipointExp::write_state` after a passed check |
| Progress params | `.param` | — (text) | `Logging::progress_save` |
| Proof points / cert | `.proof.<i>`, `.cert`, `.pack` | 6, 3, 4 | the proof tasks (`proof-system.md`) |
| Per-base children | e.g. `.div` + `add_child` names | as parent | `-divides`, `-order`, `*Generic` factor walks |

Add a state: subclass `TaskState` (framework `state-serialization.md` has the record rules), take the next free TYPE ≥ 12, update the `prst.cpp:54-65` comment and the table in §1.
