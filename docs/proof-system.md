# Proof system — deep dive

The `-proof` pipeline lets a primality test emit a small certificate that any third party can re-verify in `O(log n)` work instead of redoing the full `O(n)` exponentiation. PRST implements this in four phases: a **`SAVE`** run that records intermediate exponentiation points to disk, a **`BUILD`** run that compresses those points into a single certificate file, a **`CERT`** run that re-verifies a certificate end-to-end, and an optional **roots-of-unity** check that hardens the certificate against a specific class of forgery.

The wiring doesn't follow the normal `Run` pattern. Three things to keep in mind before reading any other code in this folder:

1. **`Proof` contains another `Run`.** Every other `Run` subclass — `Fermat`, `Pocklington`, `Morrison`, `Order`, the `*Generic` variants — is a leaf: it constructs `InputTask`s (`MultipointExp`, `LucasVMul`, …) and runs them directly. `Proof` is the only `Run` that owns a `unique_ptr<Fermat>` (`proof.h:137`) and delegates into it for SAVE/BUILD. `Pocklington` *inherits* from `Fermat` — that's normal C++; *containing* a Fermat at runtime, with ownership transferred out from under `Run::create` after the fact (`prst.cpp:360`), is not.

2. **Control flow is bidirectional.** Every other test goes top-to-bottom: `prst.cpp:419 → run->run() → task->run() → done`. SAVE/BUILD goes `prst.cpp:419 → Proof::run → Fermat::run → Proof::run(other-overload)` — Fermat calls back up into Proof at `fermat.cpp:442` to trigger the post-Fermat ProofSave/ProofBuild compression. That's why `Proof` has *two* public methods named `run` with different signatures (`proof.h:85` is the `Run::run` override; `proof.h:94` is the callback overload). It's not sloppy naming; it's the protocol.

3. **`Run::create` returns a pre-existing `Proof` for CERT.** The proof is constructed at `prst.cpp:347-348`, *before* `Run::create` is called. For CERT mode the dispatcher at `prst.cpp:458-459` short-circuits and returns that same pointer back as the Run. Every other Run subclass is allocated *inside* `Run::create`. The early construction is necessary because the `.cert` file has to be parsed (and `_Li`, `_M`, `_r_count`, `_r_0` populated) before the dispatcher could possibly know which `MultipointExp` subclass to install as `_task` for verification.

A reader who skims this code expecting the usual `Run::create → run.run() → done` pattern will get lost three different ways. Touching anywhere in the proof code without internalizing all three of the above invites silent verification holes — the kind a normal test pass won't catch.

Source files:
- `src/proof.h`
- `src/proof.cpp`
- `src/prst.cpp:344-381` (option parsing + Proof construction + ownership transfer)
- `src/fermat.cpp:240-298, 370-398, 440-444` (Fermat-side wiring of proof points and the post-Fermat callback)
- `src/exp.h:141-225, 302+` (`MultipointExp::Point`, `_on_point` callback signature, `StrongCheckMultipointExp`)
- `framework/container.h`, `framework/container.cpp` (the optional `.pack` container)
- `framework/file.h:143-158` (`FilePacked`, the container-aware `File` subclass)

Prereqs: `task-lifecycle.md` (proof tasks are `InputTask`s with non-trivial `execute()` bodies), `logging-and-progress.md` (Proof::run drives `progress().next_stage()` between phases).

## 1. The four modes at a glance

```
Proof::NO_OP = 0     no proof system involved (the default)
Proof::SAVE  = 1     -proof save N        emit N proof points to disk during a Fermat run
Proof::BUILD = 2     -proof build N       consume points, emit a .cert file
Proof::CERT  = 3     -proof cert <file>   re-verify an existing .cert
Proof::ROOT  = 4     internal-only        force the roots-of-unity check without other ops
```

Constants live at `proof.h:16-20`. `NO_OP` is the absence-marker; the other four reach `Proof::run` through different code paths. `ROOT` is not reachable via a normal user CLI invocation (there's no `-proof root` subcommand) but it is *not* dead: the `-test` harness exercises it directly. `testing.cpp:557` constructs a `Proof(Proof::ROOT, …)` and drives it through the post-Fermat overload at `testing.cpp:560`, then checks `taskRoot()->result()`: `== 1` means the residue was a non-trivial root of unity — the roots-of-unity defense caught the forged residue and the cert is rejected — so the harness flags failure only on `result() != 1` (the `Roots of unity check failed to detect the attack` path; cf. Pitfall D in §7). It's preserved as an internal trigger for the same defense `BUILD` runs by default — see §5.

A typical multi-machine workflow:

```
machine A: prst <N> -proof save 16          → produces N.proof.{0..16}
machine A: prst <N> -proof build 16         → consumes the points, produces N.cert
machine B: prst <N> -proof cert N.cert      → re-verifies the certificate cheaply
```

The same machine can do `save` and `build` back-to-back, but the modes are distinct runs because the math each performs is different. `save` is a Fermat test plus periodic disk writes; `build` is a small product-tree reduction that doesn't need to redo the Fermat work; `cert` is a single short exponentiation.

## 2. Class hierarchy

```
Proof : public Run                                        (proof.h:13)
├── Proof::Product   : TaskState   (TYPE=3)               (proof.h:23)
├── Proof::Certificate : TaskState (TYPE=4)               (proof.h:37)
└── Proof::State     : TaskState   (TYPE=6)               (proof.h:56)

ProofSave  : public InputTask                             (proof.h:140)
ProofBuild : public InputTask                             (proof.h:159)
```

`Proof` itself owns up to four sub-tasks via `unique_ptr`s:
- `_task` — `ProofSave` for SAVE, `ProofBuild` for BUILD, a `MultipointExp` subclass (`SmoothExp`/`GerbiczCheckExp`/`SlidingWindowExp`/`LiCheckExp`) for CERT.
- `_taskA` — only used in CERT mode for Li-form exponents whose top bits exceed `_M`; runs the high-bits exponentiation before the verification proper.
- `_taskRoot` — a `CarefulExp` for the roots-of-unity check, constructed in BUILD mode (and ROOT mode) when `-RootOfUnityCheck` is on (default true).
- `_fermat` — the wrapped `Fermat` instance for SAVE and BUILD; null in CERT.

The three nested `TaskState` classes are the on-disk schema for proof artifacts. **Their `TYPE` constants (3, 4, 6) are file-format identifiers — do not reuse or renumber.** The full `TYPE` registry is in `prst.cpp:54-65` and `checkpoints.md` §1.

### `Proof::State` — checkpoint shared by `ProofSave` and `ProofBuild`

```cpp
SerializedGWNum _X;        // mid-task X register (BUILD only); cheap to checkpoint
Giant            _Y;       // the current product / certificate value being built
Giant            _exp;     // running exponent for Li-mode BUILD
vector<Giant>   _h;        // hash chain — one Giant per depth level so far
```

State is keyed by `_iteration` from the base `TaskState`. For `ProofSave` it ranges from 0 to `depth()`; for `ProofBuild` it ranges from 0 to `depth() + (security ? 1 : 0)`.

### `Proof::Product` — one node of the proof-product tree

A `Giant _X` plus the `_iteration` field reused as `depth()` (the level of the product tree). Written one per file at `*.prod.{0..depth-1}` (or one entry per slot in a `.pack`).

### `Proof::Certificate` — the final artifact written by BUILD and consumed by CERT

```cpp
Giant _X;          // certificate value (a residue mod N)
Giant _a_power;    // 0 in non-Li mode; the Li-mode exponent otherwise
Giant _a_base;     // 0 in non-Li mode; the Li-mode base otherwise
```

`_iteration` is reused as `power()` (the `M` block size). The `read()` override at `proof.h:48` handles the optional `_a_power`/`_a_base` pair so an old non-Li cert remains parseable.

### `MultipointExp::Point` — proof-point schedule (lives in exp.h, not proof.h)

```cpp
struct Point {
    int  pos;       // iteration count at which this point sits
    bool check;     // whether to invoke the on_point callback here
    bool value;     // whether to capture the full Giant (vs cheap SerializedGWNum)
};
```

`Proof::calc_points` populates `proof->_points`; the vector is then handed to the `MultipointExp`/`StrongCheckMultipointExp` constructor along with the `on_point` callback so the Fermat task knows where to stop and what to emit.

## 3. Annotated `Proof::run` — the central method

`Proof::run(GWState&, File&, File&, Logging&)` (`proof.cpp:320-430`) is the `Run::run` override that gets called from `prst.cpp:419` when `-proof cert` is the active mode (or when a `-proof save/build` Run was wrapped — see §6). The body branches three ways: SAVE delegates to Fermat and tidies up, BUILD delegates to Fermat and tidies up, CERT actually does verification work.

```cpp
void Proof::run(GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging) {
    if (op() == SAVE) {
        // Hand control to the wrapped Fermat. Fermat::run sees the non-null Proof*
        // and arranges for on_point(...) to be called at each scheduled _points[i].
        _fermat->run(gwstate, file_checkpoint, file_recoverypoint, logging, this);
        _success = _fermat->success();
        _prime   = _fermat->prime();
        _res64   = _fermat->res64();

        if (_container) _container->close();        // flush the .pack file
        if (!_keep_points)
            for (int i = 1; i < count(); i++)       // 0 and count() survive — BUILD needs them
                file_points()[i]->clear();
        return;
    }
    if (op() == BUILD) {
        _fermat->run(gwstate, file_checkpoint, file_recoverypoint, logging, this);
        if (_container) _container->close();
        if (!_keep_points) {                        // BUILD output is the .cert; everything else is scratch
            if (_container) remove(_container->filename().data());
            else {
                if (!Li()) file_points()[0]->clear();   // Li-mode keeps point 0 because it's the base
                file_points()[count()]->clear();
                for (auto& f : file_products()) f->clear();
            }
        }
        return;
    }

    // CERT mode — actually do verification work.
    MultipointExp* task = dynamic_cast<MultipointExp*>(_task.get());
    if (task == nullptr) return;                    // ROOT mode falls through here; Run::create won't take you here

    logging.info("Verifying certificate of %s, complexity = %d.\n",
                 input.display_text().data(), (int)logging.progress().cost_total());
    if (gwstate.information_only) exit(PRST_EXIT_NORMAL);
    logging.set_prefix(input.display_text() + " ");

    if (Li()) {
        // Li-mode CERT: optional high-bits pre-pass with _taskA, then the strong check or plain MultipointExp.
        LiCheckExp* taskCheck = dynamic_cast<LiCheckExp*>(_task.get());
        if (taskCheck) taskCheck->init(&input, &gwstate, &file_checkpoint, &file_recoverypoint, &logging,
                                       std::move(_r_0));
        else           task->init_giant(&input, &gwstate, &file_checkpoint, &logging, std::move(_r_0));
        if (task->state() == nullptr) {
            if (_taskA) {
                // Top bits of the exponent ran first via _taskA; feed its result into the main task.
                LiCheckExp* taskACheck = dynamic_cast<LiCheckExp*>(_taskA.get());
                if (taskACheck) {
                    Giant inv_r;
                    try { inv_r = inv(_r_count, *gwstate.N); }
                    catch (const NoInverseException&) {}
                    taskACheck->init(&input, &gwstate, nullptr, nullptr, &logging,
                                     std::move(task->X0()), std::move(_r_count), std::move(inv_r));
                } else
                    _taskA->init_giant(&input, &gwstate, nullptr, &logging,
                                       std::move(task->X0()), std::move(_r_count));
                _taskA->run();
                logging.progress().next_stage();
                task->X0() = std::move(_taskA->X0());
                task->init_state(new BaseExp::StateValue(0, std::move(*_taskA->result())));
            } else {
                task->init_state(new BaseExp::StateValue(0, std::move(_r_count)));
                task->state()->set_written();
            }
        }
    } else {
        // Non-Li (smooth) CERT: just init the underlying SmoothExp / GerbiczCheckExp.
        GerbiczCheckExp* taskCheck = dynamic_cast<GerbiczCheckExp*>(_task.get());
        if (taskCheck) taskCheck->init(&input, &gwstate, &file_checkpoint, &file_recoverypoint, &logging);
        else           task->init_smooth(&input, &gwstate, &file_checkpoint, &logging);
        if (task->state() == nullptr) {
            task->init_state(new BaseExp::StateValue(0, std::move(_r_count)));
            task->state()->set_written();
        }
    }

    task->run();                                    // the real work happens here
    logging.set_prefix("");
    logging.progress().next_stage();

    _res64 = task->result()->to_res64();
    logging.result(false, "%s certificate RES64: %s, time: %.1f s.\n",
                   input.display_text().data(), _res64.data(), logging.progress().time_total());
    logging.result_save(input.input_text() + " certificate RES64: " + _res64
                        + ", time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");

    file_checkpoint.clear();
    file_recoverypoint.clear();
}
```

### The post-Fermat overload

A second `Proof::run(GWState&, Logging&, Giant* X)` overload at `proof.cpp:432-477` is **not** an override — it's a helper called by `Fermat::run` after the main Fermat exponentiation completes (`fermat.cpp:442`). Its job is to:

1. Run the roots-of-unity check on the Fermat residue `X` (if `_taskRoot` exists and `X` is non-null).
2. If `_task` is a `ProofSave`, run it now to compress the `count` proof points down to `depth` proof products.
3. If `_task` is a `ProofBuild`, run it now to compress the `depth` products down to a single certificate.

Both branches use `dynamic_cast` to pick the right mode and emit the user-visible "raw certificate RES64" / "certificate RES64" lines.

## 4. Field & method reference

### `Proof` fields (`proof.h:118-137`)

| Field           | Type                                       | Set by                       | Used by                                                    |
|-----------------|--------------------------------------------|------------------------------|------------------------------------------------------------|
| `_op`           | `int` (NO_OP/SAVE/BUILD/CERT/ROOT)         | ctor                         | every dispatch in `run()` and post-Fermat overload         |
| `_count`        | `int` (must be a power of 2 for SAVE/BUILD)| ctor                         | `calc_points`, `init_files`, file-clear loops              |
| `_Li`           | `bool`                                     | ctor (CERT) / `calc_points`  | controls smooth-vs-Li code paths everywhere                |
| `_points`       | `vector<MultipointExp::Point>`             | `calc_points`                | passed to MultipointExp; consulted by `read_point`         |
| `_M`            | `int`                                      | `calc_points` / CERT ctor    | logged via `report_param("M", _M)`; passed to certificate  |
| `_keep_points`  | `bool`                                     | `set_keep_points` (prst.cpp:373) | end-of-run cleanup branches in `run()`                |
| `_cache_points` | `bool`                                     | `set_cache_points`           | `on_point` (skip `free_buffer` if true)                    |
| `_container`    | `unique_ptr<FileContainer>`                | prst.cpp:366                 | `init_files`, end-of-run `close`, BUILD cleanup            |
| `_file_points`  | `vector<File*>`                            | `init_files`                 | `read_point`, `on_point`, cleanup                          |
| `_file_products`| `vector<File*>`                            | `init_files`                 | `read_product`, `ProofSave::execute` writes here           |
| `_file_cert`    | `File*`                                    | `init_files`                 | `ProofBuild::execute` writes; CERT ctor reads              |
| `_r_0`, `_r_count` | `Giant`                                 | CERT ctor / `init_state`     | the two endpoints of the exponentiation being verified     |
| `_r_exp`        | `Giant*`                                   | `init_state` (Li mode)       | `ProofBuild::execute` for substring decomposition          |
| `_task`         | `unique_ptr<InputTask>`                    | ctor (varies by mode)        | dispatched via `dynamic_cast` in both `run()` overloads    |
| `_taskA`        | `unique_ptr<MultipointExp>`                | CERT ctor (Li, exp > _M)     | high-bits pre-pass in `run`                                |
| `_taskRoot`     | `unique_ptr<CarefulExp>`                   | ctor (BUILD/ROOT)            | post-Fermat roots-of-unity check                           |
| `_fermat`       | `unique_ptr<Fermat>`                       | prst.cpp:360 (move-from-run) | SAVE and BUILD `run()` delegate to `_fermat->run`          |

### `Proof` methods (`proof.h:81-115`)

- `Proof(op, count, input, options, file_cert, logging)` — `proof.cpp:14`. Validates `count` is a power of 2 for SAVE/BUILD, constructs `_task` per mode, builds the optional roots-of-unity exponent for BUILD/ROOT.
- `run(gwstate, file_checkpoint, file_recoverypoint, logging)` — the `Run::run` override (`proof.cpp:320`).
- `run(gwstate, logging, X)` — the post-Fermat helper (`proof.cpp:432`); called by `Fermat::run` to run roots-of-unity, ProofSave, or ProofBuild.
- `calc_points(iterations, smooth, input, options, logging)` — `proof.cpp:121`. Builds the binary-tree `_points` schedule (see §5 lifecycle).
- `init_files(file_point, file_product, file_cert)` — `proof.cpp:181`. Allocates `count+1` point children and `depth` product children, either as plain `File`s under the parent or as `FilePacked` entries inside `_container`.
- `init_state(task, gwstate, input, logging, a)` — `proof.cpp:207`. BUILD mode: read points 0 and `count`, validate `a^k`, populate `_r_0`/`_r_count`, install initial `BaseExp::StateValue`. SAVE mode: scan `_points` from the top, pick up the highest checkpointed point as `task`'s starting state.
- `read_point(index, state, logging)` — `proof.cpp:273` and `:283` (overload taking pre-allocated state buffers). Reads `_file_points[index]`, validates `iteration() == _points[index].pos`, throws `TaskAbortException` on mismatch.
- `read_product(index, state, logging)` — `proof.cpp:295`. Same shape as `read_point` but for `_file_products`; validates `iteration() == index`.
- `on_point(index, state, logging)` — `proof.cpp:305`. The callback Fermat invokes at each scheduled point: writes the state to `_file_points[index]`, optionally drops the in-memory buffer.
- `cost()` — `proof.cpp:479`. Estimated work for progress reporting: SAVE = `(count-1)*64*1.5`, BUILD = `(depth + sec) * 2*64*1.5`, CERT = 0 (the underlying MultipointExp's cost dominates).
- Accessors: `op()`, `count()`, `Li()`, `depth()` (= `ceil(log2(_count))`), `points()`, `M()`, `container()`, `file_points()`, `file_products()`, `file_cert()`, `r_0()`, `r_count()`, `r_exp()`, `task()`, `taskRoot()`, `fermat()`. All inline.

### `ProofSave` (`proof.h:140-157`, body at `proof.cpp:488-636`)

- `init(input, gwstate, logging, proof)` — `proof.cpp:488`. Calls `InputTask::init` with `iterations = proof->count()`, sets `_state_update_period = 0` (every iteration is a checkpoint boundary — products are expensive to recompute), stashes the `Proof*`.
- `execute()` — `proof.cpp:549`. The product-tree reduction: for `i ∈ [0, depth)`, either read `_file_products[i]` if it already exists (resume case) or build the product `D` by iterating through the relevant `2^i` proof points, computing the prime-hashed combination, and writing `_file_products[i]`. Updates the running product `Y` and hash chain `h` per level. After the loop, `Y` is the raw certificate value (`Proof::State::Y()`).
- `done()` — restores logging prefix.
- `setup()` and `release()` are empty no-ops.

### `ProofBuild` (`proof.h:159-181`, body at `proof.cpp:638-771`)

- `init(input, gwstate, logging, proof)` — `proof.cpp:638`. Iterations = `proof->depth() + (security ? 1 : 0)`. If `ProofSecuritySeed` is set, derives `_rnd_seed` from the seed string concatenated with two `uint32`s of `getHighResTimer()` — non-determinism is the point of the seed (see §5).
- `execute()` — `proof.cpp:660`. Runs the Li/non-Li product-tree reduction, then optionally applies the security exponent at iteration `t` to randomize `X`/`Y`/`D` before final certificate write. Validates the final `Y != 0`. Constructs and writes a `Proof::Certificate` to `_file_cert`.
- `security()` — non-empty `_security_seed`.
- `raw_res64()` — RES64 captured *before* the security exponent is applied; used as a logged checkpoint so a third party can re-derive the same value.

## 5. Lifecycles

### SAVE — record proof points during a Fermat run

```
prst.cpp:348:   proof = new Proof(SAVE, N, ...)               // ctor reserves _task = ProofSave
prst.cpp:350:   run = Run::create(...)                        // returns a Fermat
prst.cpp:360:   proof->fermat().reset(dynamic_cast<Fermat*>(run.release()))   // move into proof
prst.cpp:372:   proof->init_files(.proof, .prod, .cert)       // allocate child Files
prst.cpp:374:   run.reset(proof.release())                    // proof is now the Run

prst.cpp:419:   run->run(...)
                ↓
proof.cpp:322:  Proof::run sees op==SAVE, delegates to _fermat->run(..., this)
                ↓
fermat.cpp:247: proof->calc_points(iterations, smooth, ...)   // build the schedule
fermat.cpp:267: bind on_point = std::bind(&Proof::on_point, proof, _1, _2, _3)
fermat.cpp:271: _task = new MultipointExp(... proof->points(), on_point)
fermat.cpp:391: proof->on_point(0, state, logging)            // write the initial value
                ↓
                Fermat exponentiation runs. At each _points[i] with check==true,
                MultipointExp::execute fires _on_point(i, state, logging) →
                Proof::on_point writes _file_points[i] (proof.cpp:305).
                ↓
fermat.cpp:442: proof->run(gwstate, logging, _task->result())
                ↓
proof.cpp:432:  post-Fermat overload runs ProofSave (compresses points → products).
                ↓
proof.cpp:329:  back in Proof::run override: close container, clear all but points 0 and count.
```

### BUILD — compress points into a certificate

Same wiring as SAVE up through `Run::create` and `init_files`. The differences:

- `_task` is a `ProofBuild`, constructed in `proof.cpp:27` with the `ProofSecuritySeed` option string.
- The Fermat run is **shorter** because the points already exist on disk; Fermat only needs to read/validate them (`init_state` at `proof.cpp:207` validates `a^k` against point 0).
- Roots-of-unity check is constructed in the Proof ctor (`proof.cpp:82-111`) and runs in the post-Fermat overload before ProofBuild (`proof.cpp:434-443`). Failure raises `TaskAbortException`.
- ProofBuild's `execute` walks each product level, hashing the products into the running `Y` value and (Li mode) accumulating `a_power`. After all levels, an optional security step applies a random exponent derived from the seed, then writes `_file_cert`.

### CERT — re-verify a certificate

Wholly separate path: no Fermat involvement, no `_fermat`, no points or products on disk other than `.cert`.

```
prst.cpp:348:   proof = new Proof(CERT, 0, ...)
                   proof.cpp:28-81: read .cert, extract _Li, _M, _r_count, _r_0
                   construct MultipointExp subclass (Smooth/GerbiczCheck or Li/SlidingWindow/LiCheck) as _task
                   if Li exponent has > _M bits, also construct _taskA for the high bits
prst.cpp:350:   Run::create returns proof (proof.cpp:458-459 in the dispatcher)
prst.cpp:419:   run->run(...) → Proof::run override falls into the CERT branch (proof.cpp:360+)
                   if Li and _taskA: run high-bits pre-pass first
                   task->run() — the actual verification exponentiation
                   _res64 = task->result()->to_res64(); print + result_save
                   clear checkpoint/recovery files
```

The verification math: CERT recomputes `X^M mod N` (non-Li) or a Li-form variant from the certificate's `r_0` to its `r_count`, where `M` is the block size from the cert. If the result matches what the cert claims, the proof is valid; the printed RES64 is the user-visible attestation.

### `calc_points` — the binary-tree schedule

`proof.cpp:121-179`. The active branch (the commented block above it is dead code) lays out `_count + 1` points across `iterations`:

- Point 0 is at iteration 0 (the base).
- Points 1..count-1 are placed at iterations forming a binary tree: the i-th point is the iteration that, when expressed in binary against `_count`, picks the path matching `i`'s bit-pattern. Concretely, the `for (int j = _count/2; j > 0 && (i & (j*2-1)) != 0; j >>= 1)` loop walks the bits of `i` and increments `pos` by the appropriate halved interval each step. This gives the canonical product-tree layout that `ProofSave::execute` then folds bottom-up.
- Point `count` is at iteration `iterations` (the final value).

Each point's `check` flag is true when `i % points_per_check == 0` (only mark a save-checkpoint at every k-th proof point to amortize disk cost); each point's `value` flag is true when it's safe to record a full `Giant` (base-2 KBNC where `k != 0`) vs. a cheaper `SerializedGWNum`.

`_M = iterations` is set on the first iteration of the building loop and then halved repeatedly; the final `_M` is what `MultipointExp` consumes as the strong-check block size.

## 6. Common patterns

### Wrapping a `Fermat` for SAVE/BUILD (`prst.cpp:344-381`)

```
1. Allocate file_cert (its TYPE = Proof::Certificate::TYPE = 4).
2. Construct Proof. Its ctor builds _task (and optionally _taskRoot) but NOT _fermat.
3. Call Run::create — for SAVE/BUILD this returns a Fermat (CERT returns the proof itself).
4. dynamic_cast<Fermat*>(run.get()) and move into proof->fermat() — this transfers ownership.
5. Allocate file_proofpoint, file_proofproduct (with proof->fingerprint(), not the input fingerprint).
6. proof->init_files(...) — wires up child Files, optionally inside the .pack container.
7. proof->set_keep_points(proof_keep) — controls cleanup at end-of-run.
8. run.reset(proof.release()) — the proof is now the active Run, and run->run will dispatch into Proof::run.
```

The fingerprint on the point/product files comes from `proof->fingerprint()` (set in the Fermat ctor at `fermat.cpp:248` to `unique_fingerprint(input.fingerprint(), "<a>.<final-pos>")`), not the raw input fingerprint. This is so different `-fermat` bases or different `-proof save N` counts get distinct point files even for the same N.

### Optional `.pack` container (`framework/container.{h,cpp}`)

When `-proof <save|build> <N> pack [filename]` is set, all proof points and products live in a single `.pack` file instead of `count+depth+1` separate files:

- The pack format is JSON-prefix-framed: a header line `["PK",1,0]\r\n`, then a sequence of `[<size>,<stream_id>]\r\n<data>` chunks, terminated by `[0,0]\r\n`. Stream id 0 is the metadata index; ids 1+ are file streams. (`container.cpp:1056, 838, 1137, 1693+`.)
- `FilePacked` (`file.h:143-158`) is a `File` subclass that delegates reads/writes through a `FileContainer`. `Proof::init_files` allocates `FilePacked` children for entries inside the container (`proof.cpp:189, 199`) and plain `File` children for entries outside.
- Corruption is detected at open time and tolerated for SAVE (resume from a partially-written pack) but warned about for BUILD (`prst.cpp:367-368`). Per-file MD5 mismatches move the corrupt entry to `_corrupted` so subsequent reads of valid entries still work.

The `.pack` form trades random-access ergonomics for one tidy file. For long SAVE runs with `count = 1024+`, this matters; for `count = 16` it's noise.

## 7. Pitfalls

### A. Don't change a `TaskState::TYPE` constant

`Proof::Product::TYPE = 3`, `Proof::Certificate::TYPE = 4`, `Proof::State::TYPE = 6`. A user with a half-completed SAVE run holds `.proof.{i}` files identified by these bytes. Bumping any of them is a contract break — tests in flight stop resuming.

If you genuinely need a schema change, add a new TYPE constant and keep the old read path. See `checkpoints.md` §1 for the registry.

### B. CERT mode bypasses `Run::create`'s normal dispatch

`Run::create` at `prst.cpp:458-459` short-circuits: if `proof->op() == Proof::CERT`, it returns the proof directly without constructing any Fermat / Pocklington / Morrison. This is why CERT skips all the half-factored-c==1 branches; the logic is "we already know what we're verifying because the cert tells us." Don't add CERT-specific dispatch elsewhere — the cert reading is in the Proof ctor, full stop.

### C. `count` must be a power of 2 for SAVE/BUILD

`proof.cpp:16` explicitly errors out if `(count & (count-1)) != 0`. The product tree assumes a balanced binary structure; an odd or non-power-of-2 count would scatter unreachable empty subtrees through `ProofSave::execute`. CERT mode doesn't use `count` (it uses `power()` from the cert), so it's exempt.

### D. The roots-of-unity check is on by default for BUILD; turning it off weakens the certificate

`proof.cpp:82` triggers `_taskRoot` construction when `op == BUILD && (!RootOfUnityCheck || RootOfUnityCheck.value())` — i.e. unless the user explicitly passes `-RootOfUnityCheck false`. The defense detects a class of forged certificates where the prover picked a clever non-trivial root of unity for the residue. `_taskRoot` is a `CarefulExp` that exponentiates the residue by the product of small prime factors of N-1; if the result is 1, the residue was a non-trivial root and the cert is rejected (`proof.cpp:438-441`).

If you turn it off you save a few seconds on BUILD; you also lose the only check that distinguishes a real proof from one constructed by an attacker who knows the factorization of N-1. **Don't disable it without understanding what you're giving up.**

### E. The security seed is not the same as the roots-of-unity check

Two different defenses, easy to conflate:
- **Roots of unity** (`-RootOfUnityCheck`, default on for BUILD): catches a forgery that uses a cleverly-chosen residue. Defense lives in `_taskRoot` which runs **post-Fermat** in the SAVE/BUILD overload at `proof.cpp:434-443`.
- **Security seed** (`-ProofSecuritySeed <string>`, default empty): catches a hypothetical bug in PRST itself that produces a certificate the verifier would re-derive identically. Defense applies a random exponent to `X`/`Y` (and `D` in Li mode) before the final cert write at `proof.cpp:736-750`. The "raw RES64" before the random exponent is logged so a third party can re-derive both values.

Use both for production-grade attestations. Use neither only for local development.

### F. `_task` polymorphism: always `dynamic_cast` before calling subclass methods

`_task` is `unique_ptr<InputTask>`. SAVE makes it a `ProofSave`, BUILD a `ProofBuild`, CERT a `MultipointExp` subclass. Code that needs subclass-specific methods (`raw_res64()`, `state()` returning `Proof::State*`, etc.) must `dynamic_cast` first. Both `run()` overloads do this — see `proof.cpp:444, 458, 360, 371, 380, 407` — and treat a null cast as "not my mode, skip this branch."

### G. `init_state` mutates the task's state out from under it

`Proof::init_state` (`proof.cpp:207`) is called from `Fermat::run` at `fermat.cpp:375` *before* `task->run()`. It either (BUILD) reads two proof points and installs a fresh `BaseExp::StateValue` as the starting state, or (SAVE) walks the points from the top picking up the most recent valid checkpoint. If you add a new mode that wraps Fermat similarly, follow this contract: install a state before `run()` or accept that the task starts at iteration 0.

## 8. Quick reference

| You want to…                                              | Call / file                                                                |
|-----------------------------------------------------------|----------------------------------------------------------------------------|
| Find which mode is active                                 | `proof->op()` returns `SAVE` / `BUILD` / `CERT` / `ROOT`                   |
| Read a proof point by index                               | `proof->read_point(i, state, logging)` (`proof.cpp:273`)                   |
| Read a product by depth-level                             | `proof->read_product(i, state, logging)` (`proof.cpp:295`)                 |
| Hook into Fermat to write a proof point                   | Already done — `Fermat::run` binds `Proof::on_point` as `_on_point`        |
| Trigger the post-Fermat ProofSave/ProofBuild compression  | Already done — `Fermat::run` calls `proof->run(gwstate, logging, X)`       |
| Force the roots-of-unity check on a non-BUILD run         | Construct `Proof(ROOT, …)` (no CLI flag exposes this today)                |
| Construct point/product files inside a .pack              | Set `proof->container()` before `proof->init_files(...)` (prst.cpp:366)    |
| Preserve intermediate files after BUILD                   | `proof->set_keep_points(true)` (prst.cpp:373; CLI: `-proof … keep`)        |
| Cache point buffers in memory between writes              | `proof->set_cache_points(true)` (skips `free_buffer` after `on_point`)     |
| Feed a verifier with a smaller working set                | Use `-proof build … -ProofPointsPerCheck k` to mark fewer `check` points   |

## 9. Open questions / non-coverage

- **The product-tree math.** This doc describes the file flow and the class shapes; the actual `hash_giants → make_prime` Fiat–Shamir construction in `ProofSave::execute` and `ProofBuild::execute` is the cryptographic core. The framework's `framework/docs/mult_*.pdf` covers the multiplication theory; the proof-system theorems themselves trace back to LiCheck (eprint 2023/195) and the Pietrzak-style sumcheck literature. Read those before modifying the math.
- **Li-mode high-bits split.** `_taskA` is constructed in CERT mode only when the Li exponent's bitlen exceeds `_M`. The choice of `SlidingWindowExp` vs. `LiCheckExp` for the high bits mirrors the main task's strong-check vs. plain mode. Worth a paragraph in a future "exponentiation algorithms" doc — `inputnum-parsing.md` and `run-hierarchy.md` are higher priority first.
- **Container internals.** `framework/container.cpp` is ~1800 lines of streaming-JSON-framed binary IO with MD5 verification, codec fields, and corruption recovery. The on-disk schema in §6 above is the contract; the implementation is documented in the framework's `container-format.md`.
- **`-proof root` would expose `Proof::ROOT`.** There's no *user* CLI surface for it; the constant exists for symmetry with the four-mode design. The constructor branch in `proof.cpp:82` (`op == BUILD || op == ROOT`) builds `_taskRoot` for it, and `testing.cpp` is a live consumer — the `-test` harness constructs `Proof(Proof::ROOT, …)` (`testing.cpp:557`) and runs it through the post-Fermat overload (`testing.cpp:560`) to verify the roots-of-unity check rejects an attacking residue. (The same harness exercises `Proof::CERT` and the post-Fermat `run()` overload too — `testing.cpp:385, 525, 549` and `testing.cpp:522, 541, 560`.) So ROOT/CERT are not unreachable; they're just not wired into a normal user invocation. To run roots-of-unity in isolation from the public CLI, both the option parser in `prst.cpp` and `Run::create`'s dispatch table would need entries.
- **`MultipointExp::Point` vs. `_points` invariants.** The `value` flag is set based on `input.type() == KBNC && k != 0 && b == 2`. The semantic is "is it cheap to materialize the full Giant at this point?", but the exact set of inputs where `value=false` is most-cost-effective is a tuning parameter not audited here. Worth a footnote in `inputnum-parsing.md`.
