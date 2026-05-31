# Task lifecycle — deep dive

The layer between `Run`-class orchestration and GWnum's FFT multiplications. A `Task` is "a long-running bignum computation that periodically checkpoints, reports progress, and survives transient FFT-numerical errors." Understanding it is the prerequisite for almost any non-trivial fix in `src/`.

Source files:
- `framework/task.h`
- `framework/task.cpp`
- `src/exp.h` / `exp.cpp` (concrete `BaseExp`-family tasks)
- `src/lucasmul.h` / `lucasmul.cpp` (Lucas-sequence tasks)
- `src/proof.h` (`ProofSave` / `ProofBuild` tasks)

Cross-references: `logging-and-progress.md` for how Task drives Progress; `lay-of-the-land.md` §"The Task layer" for one-paragraph summary.

## 1. Class hierarchy at a glance

```
TaskState                                 (base for any serializable iteration state)
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

Task                                       (abstract; restart loop, progress hook, abort)
└── InputTask                              (adds InputNum*, _timer, error-check toggles)
    ├── BaseExp                            (abstract; X_n exponentiation primitives)
    │   ├── CarefulExp                     (single careful exponentiation)
    │   ├── MultipointExp                  (multi-point exponentiation; proof-friendly)
    │   │   ├── SmoothExp                  (smooth-exponent variant: factors into small primes)
    │   │   ├── FastExp                    (binary-method exponentiation, no proof points)
    │   │   ├── SlidingWindowExp           (sliding-window exponentiation; precomputed table)
    │   │   └── StrongCheckMultipointExp   (adds Gerbicz / Gerbicz-Li check)
    │   │       ├── GerbiczCheckExp        (classic Gerbicz error check)
    │   │       ├── LiCheckExp             (Gerbicz-Li check, paper eprint 2023/195)
    │   │       └── FastLiCheckExp         (LiCheck on top of fast/sliding-window math)
    │   └── Product                        (accumulates a product of giants)
    ├── LucasVMul                          (abstract Lucas V multiplier)
    │   ├── LucasVMulFast                  (precomputed Lucas chains over a factor list)
    │   ├── LucasUVMul                     (UV form; supports strong-check)
    │   └── LucasUVMulFast                 (UV-form variant for fast paths)
    ├── ProofSave                          (-proof save: writes proof points during a Fermat run)
    └── ProofBuild                         (-proof build: consumes points, emits certificate)
```

The `TYPE` constants are on-disk discriminators in `TaskState` subclasses — they're the file format version, do not reuse or renumber.

## 2. `TaskState` — the serializable iteration anchor

```cpp
class TaskState {
    char _type;           // TYPE constant, sets the on-disk version
    bool _written;        // whether this state has been flushed to _file already
    int  _iteration;      // current iteration count (the only field the base class persists)

    void set(int iteration);      // resets _written to false and stores _iteration
    virtual bool read(Reader&);   // base reads _iteration only
    virtual void write(Writer&);  // base writes _iteration only
    void set_written();           // mark state as flushed; Task::write_state skips re-writing
    char type();                  // returns _type (the on-disk discriminator)
    char version() { return 0; }  // unused on-disk version field; reserved for future schema bumps
    bool is_written();             // returns _written
    int  iteration();              // returns _iteration
};

template<class State>
State* read_state(File* file);   // null-safe deserializer
```

Subclasses override `read()` / `write()` and add their own fields (`_giant_value`, `_serialized_value`, `_V` / `_Vn` / `_Vn1`, etc.). `read_state<T>(file)` is the canonical way to load a checkpoint: returns a fresh `T` if the file decodes, or `nullptr`. The `_written` flag prevents re-writing a state that's already on disk in `Task::write_state()` (`task.cpp:192`).

Subclasses' `set(...)` methods take the same `int iteration` first arg plus any state-specific payload. `BaseExp::StateValue::set(iteration, GWNum&)` for example writes a `Giant`; `BaseExp::StateSerialized::set(iteration, GWNum&)` writes a `SerializedGWNum` (cheaper to checkpoint mid-task).

## 3. `Task` fields — annotated

```cpp
class Task {
    bool _error_check = false;                                  // promoted to ReliableGWArithmetic when true
    arithmetic::GWState* _gwstate = nullptr;                    // owning GWnum runtime; configured by main()
    arithmetic::GWArithmetic* _gw = nullptr;                    // active arithmetic; flips between two slots in run()
    std::unique_ptr<TaskState> _state;                          // current durable state
    std::unique_ptr<TaskState> _tmp_state;                      // scratch slot, swapped into _state to amortize allocation
    int _iterations = 0;                                        // total iteration count for this task instance
    int _state_update_period = MULS_PER_STATE_UPDATE;           // how often commit_execute checkpoints
    File* _file;                                                // checkpoint File (may be null for transient tasks)
    Logging* _logging;                                          // Logging or SubLogging instance for progress + reports
    std::chrono::system_clock::time_point _last_write;          // for DISK_WRITE_TIME pacing
    std::chrono::system_clock::time_point _last_progress;       // for PROGRESS_TIME pacing
    static bool _abort_flag;                                    // process-wide cooperative cancel
    int _restart_count = 0;                                     // restarts within current execute() pass
    int _restart_op = 0;                                        // op index from which to resume after a reliable restart
    double _op_count = 0;                                       // accumulated op count across restarts
    double _op_base = 0;                                        // _gwstate->ops() snapshot at start of current pass
};
```

Public surface:

```cpp
static void abort();              // signal handler calls this from prst.cpp:38
static void abort_reset();        // clear the static _abort_flag (used by tests / batch mode)
static bool abort_flag();         // read the static _abort_flag
virtual void run();               // the main entry point
arithmetic::GWArithmetic& gw();   // active arithmetic (after init); used by setup()/execute() to issue ops
TaskState* state();               // current durable state; null before first commit
int iterations();                 // total iteration count for this task instance (set by init)
int state_update_period();        // _state_update_period; commit_execute uses it to gate checkpoints
bool is_last(int iteration);      // true when next commit will checkpoint
virtual double progress();        // 0..1; default = state()->iteration() / iterations()
int ops();                        // accumulated reliable op count
```

### `InputTask` adds

```cpp
class InputTask : public Task {
    InputNum* _input = nullptr;         // candidate being tested; needed for FFT re-setup on restart
    double _timer = 0;                  // wall-clock anchor; getHighResTimer() at init, elapsed at done()
    bool _error_check_near = true;      // auto-enable error check when near FFT limit
    bool _error_check_forced = false;   // unconditionally enable error check

    void set_error_check(bool near, bool check);   // toggle the two flags above; recomputes _error_check
    double timer();                     // !! see Pitfall G — only meaningful AFTER done()
};
```

Most concrete tasks derive from `InputTask`. The `_timer` field is the wall-clock measurement of how long this single task ran.

## 4. The cadence knobs

Three static knobs on `Task` control pacing (`task.cpp:27-29`):

| Constant                   | Default | What it gates                                                                                       |
|----------------------------|---------|-----------------------------------------------------------------------------------------------------|
| `MULS_PER_STATE_UPDATE`    | 20000   | Default `_state_update_period`. Subclasses may shrink it (e.g. `MultipointExp` for high-base smooth exponents at `exp.cpp:135`). |
| `DISK_WRITE_TIME`          | 300 s   | Wall-clock interval between disk checkpoints in `on_state()`.                                       |
| `PROGRESS_TIME`            | 10 s    | Wall-clock interval between `report_progress()` calls in `on_state()`.                              |

`MULS_PER_STATE_UPDATE` is iteration-based; the others are time-based. They run independently. Both `DISK_WRITE_TIME` and `PROGRESS_TIME` use `std::chrono::system_clock` (not `getHighResTimer`).

## 5. Annotated `Task::run()`

The retry-loop is the hardest part to read at a glance. From `task.cpp:45-151`, with annotations:

```cpp
void Task::run() {
    ReliableGWArithmetic* reliable;
    std::unique_ptr<arithmetic::GWArithmetic> arithmetics[2];   // [0]=setup pass, [1]=execute pass

    int setup_restart_count = 0;
    _restart_count = 0;
    _restart_op = 0;
    int* restart_count[2] { &setup_restart_count, &_restart_count };

    while (true) {                                              // OUTER: per-FFT-size attempt
        _op_base = _gwstate->ops();
        int i;
        for (i = 0; i < 2; ) {                                  // INNER: 0=setup, 1=execute
            if (!arithmetics[i])
                arithmetics[i].reset(_error_check ? new ReliableGWArithmetic(*_gwstate)
                                                  : new GWArithmetic(*_gwstate));
            _gw = arithmetics[i].get();
            reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);

            try {
                if (i == 0) setup();
                if (i == 1) execute();
                *restart_count[i] = 0;                          // pass succeeded; clear its counter
                i++;
                continue;
            }
            catch (const TaskRestartException&) {
                if (!reliable) {                                // first restart: promote to reliable
                    release(); _gw = nullptr;
                    _op_count = _gwstate->ops() - _op_base;
                    _op_base = _gwstate->ops();
                    _error_check = true;                        // upgrade flag, drop arithmetics
                    arithmetics[0].reset(); arithmetics[1].reset();
                }
                // else: fall through to restart logic below
            }
            catch (const ArithmeticException& e) {
                _logging->warning("ArithmeticException: %s\n", e.what());
            }
            catch (...) {
                release(); throw;                                // unrecoverable, bubble up
            }

            _logging->warning("Arithmetic error, restarting at %.1f%%.\n", 100.0*progress());

            if (reliable && reliable->restart_flag() && !reliable->failure_flag()) {
                // reliable detected a suspicious op; restart from _restart_op
                std::string opstr = "pass " + std::to_string(i) + " suspicious ops:";
                for (auto it = reliable->suspect_ops().begin(); it != reliable->suspect_ops().end(); it++)
                    if (*it >= (i == 1 ? _restart_op : 0))
                        opstr += " " + std::to_string(*it);
                opstr += ".\n";
                _logging->debug(opstr.data());
                reliable->restart(i == 1 ? _restart_op : 0);
                i = 0;                                          // re-run setup before re-execute
                continue;
            }
            if ((!reliable || !reliable->failure_flag()) && *restart_count[i] < 2) {
                if (reliable) reliable->restart(i == 1 ? _restart_op : 0);
                (*restart_count[i])++;
                i = 0;
                continue;
            }
            break;                                              // restart logic exhausted; bump FFT
        }
        if (_gw != nullptr) _op_count = _gwstate->ops() - _op_base;
        if (i == 2) break;                                      // both passes succeeded; we're done

        // ESCAPE: bump FFT and retry from scratch
        if (_gw != nullptr) release();
        _gw = nullptr;
        if (_gwstate->next_fft_count >= 5) {
            _logging->error("Maximum FFT increment reached.\n");
            throw TaskAbortException();                          // give up
        }
        if (*restart_count[i] >= 5) {
            _logging->error("Too many restarts, aborting.\n");
            throw TaskAbortException();
        }
        (*restart_count[i])++;
        _gwstate->next_fft_count++;
        _logging->report_param("next_fft", _gwstate->next_fft_count);
        reinit_gwstate();                                        // virtual; InputTask reconfigures FFT
        _restart_op = 0;
        arithmetics[0].reset(); arithmetics[1].reset();
    }
    _tmp_state.reset();
    if (_gw != nullptr) release();
    _gw = nullptr;
}
```

Mental model:

- **Two passes.** `setup()` runs once before `execute()`. They use separate `GWArithmetic` instances (slot 0 vs slot 1). Some subclasses do non-trivial work in `setup()` (e.g. `MultipointExp::setup` builds sliding-window tables); others leave it nearly empty.
- **Three escape routes.** (a) Both passes succeed → break the outer loop. (b) Reliable detected a suspect op → re-restart from `_restart_op` without bumping FFT. (c) Restart budget exceeded → bump `next_fft_count`, call `reinit_gwstate()`, retry. After 5 FFT increments or 5 same-pass restarts → `TaskAbortException`.
- **Promotion to reliable.** First `TaskRestartException` upgrades the task: `_error_check = true`, drop both arithmetic slots, the next iteration constructs `ReliableGWArithmetic`. After that the suspect-op restart path is available.

## 6. Annotated `Task::on_state()`

The periodic callback every concrete task hits via `commit_execute<TState>`. From `task.cpp:163-190`:

```cpp
void Task::on_state() {
    bool state_save_flag = _logging->state_save_flag();

    // 1. Disk-save tick: every DISK_WRITE_TIME seconds, on abort, or on logging request.
    if (now - _last_write >= DISK_WRITE_TIME || abort_flag() || state_save_flag) {
        _logging->debug("saving state to disk.\n");          // visible only at LEVEL_DEBUG
        _logging->progress().update(progress(), ops());     // anchors timer / accrues elapsed
        _logging->state_save();                              // virtual; default no-op, BOINC overrides
        write_state();                                       // _file->write(*_state) if file present
        _last_write = now;                                   // reset DISK_WRITE_TIME pacer
        _logging->progress().update(progress(), ops());     // re-anchor after I/O wall-clock
        _logging->progress_save();                           // .param file write
    }

    // 2. Honor cooperative abort BEFORE potentially running progress logic.
    if (abort_flag()) throw TaskAbortException();            // SIGTERM/SIGINT path; caught in run() callers

    // 3. Progress tick: every PROGRESS_TIME seconds.
    if (now - _last_progress >= PROGRESS_TIME) {
        _logging->progress().update(progress(), ops());     // accrue elapsed before formatting
        _logging->report_progress();                         // emits "X.X% done…" line
        _last_progress = now;                                // reset PROGRESS_TIME pacer
    }

    // 4. If we're running reliable, snapshot the latest known-good op for restart.
    if (_error_check && _gw != nullptr) {
        ReliableGWArithmetic* reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
        _restart_op = reliable->op();
    }

    // 5. Always: heartbeat (BOINC).
    _logging->heartbeat();
}
```

Cross-ref: this is the place inside the framework that drives the **recurring/periodic** `Progress::update()` ticks from a Task's perspective. It is not the only call site, though: `Task::init` calls `progress().update(0, ops())` to anchor the timer at task start (`task.cpp:40`), and `InputTask::done` calls `progress().update(1, ops())` to mark the stage 100% complete at task end (`task.cpp:235`). Every long task automatically generates the periodic ticks that populate `_time_total` on its `Logging`'s `Progress` (see `logging-and-progress.md` §7).

## 7. State machinery

The pattern is:

```cpp
// In subclass execute() loop:
for (i = state ? state->iteration() : 0; i < iterations(); i++) {
    /* one iteration of math: e.g. X = X*X */
    commit_execute<MyState>(i + 1, /* ctor args for MyState::set */);
}
```

What `commit_execute<T>` does (`task.h:96`):

1. Checks the cadence condition: `iteration - last_state_iter >= state_update_period || iteration == iterations() || abort_flag() || state_save_flag()`.
2. If true: call `check()` (which throws `TaskRestartException` if reliable detected a problem), then `set_state<T>(iteration, args...)`.
3. If false: skip — let the math loop continue without snapshotting.

What `set_state<T>` does (`task.h:105`):

1. If `_tmp_state` is null or wrong type, allocate a fresh `T` into `_tmp_state`.
2. Call `_tmp_state->set(iteration, args...)`.
3. **Swap** `_tmp_state` and `_state`. Now the new state is durable, the old one is the scratch slot for next iteration. This is why two slots exist — to amortize allocation across iterations.
4. Call `on_state()` — which might checkpoint to disk, report progress, throw on abort, or do nothing.

`reset_state<T>()` (`task.h:113-118`) would drop any existing state, install a fresh `T()`, and call `on_state()` — but it has **no live callers** anywhere in the tree, so treat it as effectively dead code rather than a recommended pattern. When a subclass wants to start clean (e.g., `MultipointExp::execute` at `exp.cpp:202-206` when there's no checkpoint and no `_smooth`), the real pattern is to construct the state directly and hand it to `init_state` — `tmp_state = new StateSerialized()` followed by `init_state(tmp_state)` — not `reset_state`.

`is_last(iteration)` (`task.h:82`) is the same predicate as `commit_execute`'s gate. Subclasses use it to decide whether the next iteration is a snapshotting boundary — useful when the math has work that should only happen at boundaries (e.g. switching from `SerializedGWNum` to a final `Giant`).

`commit_setup()` (`task.cpp:203`) is a different call, made from inside `setup()`: it runs `check()` and resets reliable state if needed. Subclasses that perform work in `setup()` should call it at the end of that work.

## 8. The error-correction subsystem

Two exception types drive recovery:
- `TaskRestartException` — "this pass needs to redo from a known-good op." Caught in `run()`'s inner loop.
- `TaskAbortException` — "abandon this task entirely." Caught by callers (e.g., `prst.cpp:423`).

`ReliableGWArithmetic` is GWnum's verifier-wrapping arithmetic. It exposes:
- `restart_flag()` — "I detected something suspicious; please restart."
- `failure_flag()` — "the suspicious op was confirmed bad; restart from a backup."
- `suspect_ops()` — the op indices that look wrong.
- `restart(op)` — rewind internal state to operation `op`.

`Task::check()` (`task.cpp:153`) is called from inside `commit_execute` and `commit_setup`. If reliable's `restart_flag` or `failure_flag` is set, it throws `TaskRestartException`. The exception is then caught in `run()`, which decides whether to (a) re-run from `_restart_op` (still inside the same FFT size) or (b) bump FFT.

`_restart_op` is updated inside `on_state()` (`task.cpp:184-188`): when running reliable, the most recent reliable op count is snapshotted there. After a restart, that's where execution resumes.

`InputTask::reinit_gwstate()` (`task.cpp:221`) is what `Task::run` calls when bumping FFT size: tears down `_gwstate`, calls `_input->setup(*_gwstate)` again, warns the user (preserving any prefix), recomputes `_error_check_near` against the new FFT.

## 9. `InputTask::init` and `done`

```cpp
void InputTask::init(InputNum* input, GWState* gwstate, File* file, TaskState* state, Logging* logging, int iterations) {
    Task::init(gwstate, file, state, logging, iterations);
    _input = input;                                              // remember for FFT re-setup on restart
    _timer = getHighResTimer();                                  // start wall-clock anchor
    _error_check = _error_check_near                             // either: auto-on near FFT limit
        ? gwnear_fft_limit(gwstate->gwdata(), 1) == TRUE         // or: respect the forced flag
        : _error_check_forced;
}

void InputTask::done() {
    _timer = (getHighResTimer() - _timer) / getHighResTimerFrequency();   // convert to elapsed seconds
    _logging->progress().update(1, ops());                       // mark stage 100% complete; flushes elapsed up the parent chain
    _logging->debug("task time: %.3f s, ops: %d, time per op: %.3f ms.\n",  // debug-only summary
                    _timer, ops(), _timer*1000/ops());
}
```

Key consequence: `task->timer()` is **the start anchor** during `run()`, and **the elapsed seconds** only after `done()` returns. Pitfall G below.

`_error_check` is decided at init based on whether the FFT is close to its safety limit (`_error_check_near` mode) OR forced on (`_error_check_forced`). Subclasses can change this with `set_error_check(near, check)` before running.

## 10. Concrete subclasses — one-paragraph tour

- **`BaseExp`** (`exp.h:11`). Abstract base for X^k mod N exponentiation. Owns `_exp` (the exponent), `_tail`, `_X0` / `_x0` (starting value), `_smooth` flag. Provides three nested state classes: `State` (abstract), `StateSerialized` (cheap mid-task representation), `StateValue` (final-iteration `Giant`). The override of `commit_execute` at `exp.h:67-74` switches automatically between the two state types based on whether `iteration == iterations()`.

- **`CarefulExp`** (`exp.h:93`). Single careful exponentiation; no error-check restart path. Used for short side-tasks (Pocklington's per-factor checks, Order's tail).

- **`MultipointExp`** (`exp.h:138`). Exponentiation that hits a list of `_points` (intermediate iteration counts) where checkpoints are forced. Used by the proof system (each point is a proof point), and by Fermat tests with checkpoints. Shrinks `_state_update_period` for high-base smooth exponents (`exp.cpp:135`).

- **`StrongCheckMultipointExp` / `GerbiczCheckExp` / `LiCheckExp` / `FastLiCheckExp`** (`exp.h:302+`). Layer Gerbicz / Gerbicz-Li strong error checks on top of MultipointExp. They use `StrongCheckState` to persist the check-side intermediate. The math is documented in the framework's `docs/mult_*.pdf`.

- **`Product`** (`exp.h:545`). Accumulates a product of giants. Used in `MorrisonGeneric` / `PocklingtonGeneric` to compute the product of cofactor pieces (`morrison.cpp:531`, `pocklington.cpp:377`).

- **`LucasVMul` / `LucasVMulFast`** (`lucasmul.h:9, 29`). Lucas V-sequence multiplier. State holds a single `Giant V` plus a parity flag and an index. `LucasVMulFast` exposes `mul_giant` / `mul_prime` to register the factor list before running.

- **`LucasUVMul` / `LucasUVMulFast`** (`lucasmul.h:104, 206`). UV-form Lucas multiplier, with an optional strong-check via `StrongCheckState`. Used by `Morrison`'s main test path.

- **`ProofSave` / `ProofBuild`** (`proof.h:140, 159`). Proof-system tasks: `ProofSave` writes proof points during a Fermat run; `ProofBuild` consumes them on the verifier side. Both layer their own state types over `InputTask`.

## 11. Common patterns

### Pattern: a Run-class running a single task

```cpp
// In Foo::run()
File* checkpoint = file_checkpoint.add_child(...);
_task->init(&input, &gwstate, checkpoint, &logging, ...);
_task->run();
Giant* result = _task->result();
if (result == nullptr) { /* aborted or failed */ }
// ... use result ...
_task.reset();
checkpoint->clear();
```

After `run()` returns normally, `result()` is non-null and `task->timer()` is the elapsed seconds.

### Pattern: resume from checkpoint

`init()` is responsible for picking up state. The convention:

```cpp
void Foo::init(...) {
    InputTask::init(input, gwstate, file, /*state*/ nullptr, logging, iterations);
    State* state = read_state<State>(file);                  // null-safe load
    if (state != nullptr) init_state(state);                 // restore _state, sync progress
}
```

If the file is fresh, `read_state` returns null and the task starts at iteration 0. If it has saved data, the task picks up at `state->iteration()`.

### Pattern: sub-task under a `SubLogging`

See `logging-and-progress.md` §9 for the full pattern. Briefly: register an outer stage, hand the inner task a `SubLogging`, run it, then call the **outer** `progress().update(progress_stage(), 0)` to credit elapsed time durably to outer's `_time_total`.

### Pattern: chained tasks (Fermat → factor checks)

`Fermat::run` runs the main exponentiation, then `Pocklington::run` (or its loop) creates per-factor `CarefulExp` tasks one at a time, each fed by the running result of the previous. Each sub-task has its own `init` / `run` / `result` lifecycle; checkpoints are scoped via `file.add_child(...)`.

## 12. Pitfalls

### A. Subclass `setup()` without `commit_setup()` at the end
`commit_setup()` runs `check()` and resets reliable state. If your `setup()` does GW operations and skips `commit_setup`, a suspicious-op detection during setup never escalates into a `TaskRestartException`, and `_restart_op` semantics get confused. Rule: any `setup()` that touches `_gw` must end with `commit_setup()`. See `MultipointExp::setup` at `exp.cpp:171-172` for the canonical pattern.

### B. Mutating `_state` directly
`_state` is a `unique_ptr<TaskState>`. Don't reassign it from a subclass — go through `set_state<T>` / `commit_execute<T>` (or `reset_state<T>`, which routes through `on_state()` correctly but has no live callers — see §7). Direct assignment bypasses `on_state()`, so no progress tick, no checkpoint, no abort handling, and the `_tmp_state` swap protocol breaks.

### C. Reading `task->result()` before completion
`BaseExp::result()` returns `nullptr` if `state()->iteration() != iterations()` or if the state isn't a `StateValue` (`exp.h:66`). `LucasVMulFast::result()` and `LucasVMul::result()` similarly guard. Don't dereference before `run()` returns.

### D. Assuming `_state_update_period == MULS_PER_STATE_UPDATE`
`MultipointExp::init` (`exp.cpp:133-135`) shrinks `_state_update_period` for smooth exponents with base > 2. If you write code that's tuned to the static default (e.g. estimating disk bandwidth from a known checkpoint cadence), it'll be wrong for those cases. Read `task->state_update_period()` if you need the actual value.

### E. Calling `Task::abort()` from a non-signal-handler thread
Currently fine — it just sets `bool _abort_flag`. But it isn't documented as MT-safe; if anyone ever runs Tasks in parallel above GWnum, this needs `std::atomic<bool>` and acquire/release semantics. Treat `abort()` as a process-wide one-shot.

### F. Subclass swallowing `TaskRestartException` / `TaskAbortException`
Both are caught by `Task::run()` and are part of its protocol. If a subclass `setup()`/`execute()` catches them locally and converts them into something else, the restart machinery never engages. Rule: don't catch `TaskRestartException` or `TaskAbortException` inside `setup` / `execute`; let them propagate.

### G. Reading `InputTask::timer()` mid-run
During `run()`, `_timer` holds the start `getHighResTimer()` reading. Only after `done()` does it become elapsed seconds. If you log or push `task->timer()` from inside the run loop, you'll get garbage (a ~13-digit timer count). Read it only after the task completes — typically right after the next `task->run()` call returns or in a result-gathering pass.

### H. Forgetting `add_child` for sub-task files
A sub-task that shares the parent's checkpoint file will overwrite the parent's state. Always allocate a child file via `file.add_child(name, fingerprint)` for any sub-task. See `morrison.cpp:224` for the canonical uniqueness pattern (`checkpoint = file_checkpoint.add_child(sP, File::unique_fingerprint(file_checkpoint.fingerprint(), sP));`), with the recovery-point analogue at `:229`; `File::unique_fingerprint` is also used at `:363` and `:516`.

## 13. Quick reference

| You want to…                                                        | Call                                                                                  |
|---------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| Define a new persistent state                                       | Subclass `TaskState`, pick a unique `TYPE`, override `read`/`write`                   |
| Define a new task                                                   | Subclass `InputTask`, override `setup` + `execute` + `release` (+ `reinit_gwstate` if non-trivial) |
| Advance state inside an iteration loop                              | `commit_execute<MyState>(iter, args...)` after each math step                         |
| Force a checkpoint to disk                                          | `commit_execute<MyState>(iterations(), args...)` (the `iter == iterations()` path)    |
| Start state clean (no checkpoint)                                   | `new MyState()` + `init_state(state)` (see `exp.cpp:202-206`); `reset_state<MyState>()` exists at `task.h:113-118` but has no live callers |
| Read state from a checkpoint file                                   | `read_state<MyState>(_file)` in `init`, then `init_state(state)`                      |
| Signal a restart from a known-good op (reliable mode)               | Throw `TaskRestartException` (or rely on `check()` to throw via reliable's flags)     |
| Cooperatively cancel from outside                                   | `Task::abort()` (sets the static bool; tasks honor it at next `on_state`)             |
| Capture a task's wall-clock duration                                | Read `((InputTask*)task)->timer()` **after** `task->run()` returns                    |
| Allocate a unique sub-task checkpoint file                          | `parent_file.add_child(name, File::unique_fingerprint(parent_fp, name))`              |
| Force error-check on for the whole run                              | `task->set_error_check(/*near*/false, /*check*/true)` before `task->init`             |
| Bump checkpoint frequency for a single task                         | Set `_state_update_period` after calling base `init`                                  |

## 14. Open questions / non-coverage

- **The math.** Gerbicz-Li's algorithm, Lucas chain construction, sliding-window exponentiation, Proth/Pocklington/Morrison theorem details. See `framework/docs/mult_en_20230925.pdf` for the multiplication theory; the test theorems themselves are in the upstream README and the original papers.
- **GWnum internals.** FFT selection, thread-pool layout, instruction-set dispatch. Owned by `framework/gwnum/` (third-party).
- **`ReliableGWArithmetic` internals.** How GWnum decides an op is suspect, what `failure_flag` means concretely. Treat it as a black box for now — its public flags are the contract.
- **State serialization format.** `Giant`/`SerializedGWNum`/etc. all have custom `read`/`write` that go through `Reader`/`Writer` in `framework/file.cpp`. Worth its own short doc if we ever need to bump a `TaskState::TYPE` constant or extend the on-disk schema.
- **The `_smooth` exponentiation path.** `BaseExp::_smooth` toggles a different math path (used when the exponent factors smoothly into small primes). Worth a footnote in any future "exponentiation algorithms" doc.
