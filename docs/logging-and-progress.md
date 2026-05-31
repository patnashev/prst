# Logging & Progress — deep dive

One of the subtlest subsystems in the codebase. This doc is the reference: what each field means, what each method does, the invariants you must preserve, and the patterns/pitfalls that come from misreading them.

Source files:
- `framework/logging.h`
- `framework/logging.cpp`
- `framework/task.cpp` (only the call sites that drive Progress)

## 1. The three classes at a glance

```
Progress      ← pure data + arithmetic; no I/O. Owns staged costs and time bookkeeping.
Logging       ← owns one Progress + result/factor/log files + the progress checkpoint file.
SubLogging    ← inner-scope Logging that forwards most calls to a parent Logging.
```

Each `Logging` has exactly one `Progress` (`Logging::_progress`). `SubLogging` is a `Logging` whose `_progress._parent` is set to the parent `Logging`'s `Progress`. There is no shared global state; the only coupling between sub and parent is the parent pointer chain inside `Progress`.

## 2. Progress: fields and what they mean

```cpp
class Progress {
    std::vector<double> _costs;     // per-stage expected work (registered via add_stage)
    double _total_cost = 0;         // sum of _costs
    int    _cur_stage = 0;          // index into _costs
    double _cur_progress = 0;       // 0..1, fractional progress through _costs[_cur_stage]
    double _time_total = 0;         // DURABLE elapsed seconds; only update()/time_init() change it
    double _time_stage = 0;         // TRANSIENT; reset to 0 by every update() on this Progress
    double _time_op = 0;            // most recent estimated seconds-per-operation
    double _timer = 0;              // last getHighResTimer() reading; 0 ⇒ uninitialized
    int    _op_count = 0;           // most recent ops() snapshot
    Progress* _parent = nullptr;    // set by set_parent(); SubLogging wires this in its ctor
    std::map<std::string, std::string> _params;   // name=value pairs persisted to .param file
};
```

**The two time fields are the trap.** They serve different roles:

| Field         | Lifecycle                                                                  | Who increments it                                                                     |
|---------------|----------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| `_time_total` | Persists across `update()` calls and across runs (saved to `.param`).      | `update()` on **this** Progress; `time_init()` on resume; nothing else.               |
| `_time_stage` | Transient — wiped to 0 by every `update()` on **this** Progress.           | `update()` on a **descendant** Progress (the parent walk adds elapsed up the chain).  |

`time_total()` returns `_time_total + _time_stage`, so a *correct* report at any moment is the sum of "what's been credited durably" + "what descendants have added since this Progress's own last update."

## 3. Progress: the methods

```cpp
void reset();                                  // zero all bookkeeping
void add_stage(double cost);                   // append to _costs, increase _total_cost
void next_stage();                             // _cur_stage++; then update(0, 0)
void skip_stage();                             // _cur_stage++ only — does NOT call update()
void update(double progress, int op_count);    // see §4
void time_init(double elapsed);                // _time_total = elapsed; _timer = now
void set_parent(Progress* parent);             // wires _parent for the update() walk; SubLogging ctor uses it

double progress_stage();                       // _cur_progress
double progress_total();                       // weighted progress across all stages, 0..1
double cost_total();                           // _total_cost
double time_total();                           // _time_total + _time_stage   ← composed
double time_op();                              // _time_op
int    op_count();                             // _op_count (last snapshot passed to update())
int    num_stages();                           // _costs.size()
int    cur_stage();                            // _cur_stage

std::string& param(const std::string& name);   // mutable; auto-creates
int    param_int(const std::string& name);     // strtol of the string param ("" → 0)
double param_double(const std::string& name);  // strtod of the string param ("" → 0.0)
std::map<std::string, std::string>& params();  // raw map; used by progress_save / file_progress
```

The correct approach is to call `update()` on the outer Progress with `(progress_stage(), 0)`; that captures wall-clock into `_time_total` via the existing `_timer` anchor. See §6 and §10.

## 4. `update()` — the heart of the mechanism

This is the only function you really need to internalize. From `framework/logging.cpp:10-31`:

```cpp
void Progress::update(double progress, int op_count) {
    _cur_progress = progress;
    if (_timer == 0)
        _timer = getHighResTimer();
    double elapsed = (getHighResTimer() - _timer) / getHighResTimerFrequency();
    _timer = getHighResTimer();
    _time_total += elapsed;        // 1. credit wall-clock since last update on THIS progress
    _time_stage = 0;               // 2. wipe transient field
    if (_op_count < op_count)
        _time_op = elapsed / (op_count - _op_count);
    _op_count = op_count;

    Progress* cur = this;
    Progress* parent = _parent;
    while (parent != nullptr) {    // 3. walk UP, poking each ancestor
        parent->_cur_progress = cur->progress_total();
        parent->_time_stage  += elapsed;
        cur = parent;
        parent = cur->_parent;
    }
}
```

What this does, in three steps:

1. **Credit elapsed time durably on this Progress.** `_time_total += (now − _timer)`. `_timer` is the wall-clock anchor; once initialized it's monotonic.
2. **Reset `_time_stage` on this Progress to 0.** Important: the design rationale is that `_time_stage` is "what descendants have added since *my* last update." Once *I* have updated, that contribution has been observed, and starts over.
3. **Walk the parent chain** updating each ancestor: their `_cur_progress` is set to this descendant's overall progress, and the same `elapsed` is added to each ancestor's `_time_stage`. (Note: NOT to `_time_total`. That's deliberate — only an ancestor's own `update()` can credit `_time_total`.)

Implications worth pinning:

- **Calling `update()` on a SubLogging's Progress ≠ calling `update()` on the parent.** The sub's update increments parent's `_time_stage` only. The parent's `_time_total` is not touched until *the parent's own* `update()` fires.
- **The first call ever to `update()` returns `elapsed = 0`** (because `_timer == 0` initially gets initialized to `now`). This is why `Logging::file_progress()` calls `time_init()` at startup — it pre-arms `_timer` so the first real `update()` is meaningful.
- **`_time_stage` is not additive across siblings.** If two children both call `update()`, the parent's `_time_stage` gets `elapsed_1 + elapsed_2`. But the next call to the parent's own `update()` will reset it to 0, capturing the cumulative wall-clock into `_time_total` instead.

## 5. Logging: fields and methods

```cpp
class Logging {
    int     _level;                      // gate for debug/info/warning/error/result
    Progress _progress;                  // the single owned Progress; exposed via progress()
    File*   _file_progress = nullptr;    // .param checkpoint
    std::string _file_result = "result.txt";       // result lines persisted here via result_save()
    std::string _file_factor = "factors.txt";      // factors emitted by report_factor() appended here
    std::string _file_log;               // optional log file
    std::string _prefix;                 // printed before every report() line
    bool    _print_prefix = true;        // controls whether _prefix is printed next; toggled by report() based on prior message
    bool    _overwrite_line = false;     // for `\r`-style progress overwrite
};
```

Levels (sorted ascending):

```
LEVEL_DEBUG   = 0
LEVEL_INFO    = 1
LEVEL_PROGRESS= 2     ← report_progress() emits at this level (overwriting line)
LEVEL_WARNING = 3
LEVEL_ERROR   = 4
LEVEL_RESULT  = 5     ← result() emits at this level
```

`report()` is the single output funnel. `info/warning/error/debug/result` all build a buffer with `vsnprintf`, gate against `_level`, and call `report()`. Two side effects in `report()`:
- Optional progress-line overwrite: if the previous output was a progress line (`_overwrite_line == true`), prefix the new line with a CR + spaces to overwrite it.
- Optional log-file mirror: if `_file_log` is set and the level is not `LEVEL_PROGRESS`, append to that file.

`warning()` and `error()` auto-persist: each calls `result_save(buf)` immediately after `report()` (`logging.cpp:72-73` and `logging.cpp:85-86`). `result()` does **not** auto-persist — its body (`logging.cpp:89-99`) only calls `report()`. Result lines land in `_file_result` only because callers invoke `logging.result_save(...)` separately, often with a slightly different string (as the §11 table shows, and as at `fermat.cpp:420-421`, where `result()` prints `time_total()` as a float and the following `result_save()` writes the integer-second form).

`report_progress()` formats one of the two human-readable progress strings (with or without "stage / total" depending on `num_stages() > 1`) and emits at `LEVEL_PROGRESS`. It is called from `Task::on_state()` every `Task::PROGRESS_TIME` seconds (default 10).

`progress_save()` writes `_progress._params` (including the `time_total` parameter) to the `.param` file. `Task::on_state()` triggers this every `Task::DISK_WRITE_TIME` seconds (default 300) and also on abort. Empty values are skipped.

`file_progress()` reads the `.param` file at startup, repopulates `_params`, and calls `_progress.time_init(_params["time_total"])`. This is how time accounting survives restarts.

## 6. SubLogging: what it overrides

```cpp
class SubLogging : public Logging {
    SubLogging(Logging& parent, int level = LEVEL_WARNING)
        : Logging(level), _parent(parent) {
        progress().set_parent(&parent.progress());     // <-- the parent walk hookup
    }

    // forward-to-parent overrides:
    report(message, level)            → _parent.report(prefix+message, level)   // (prefix added once)
    report_param(name, value)         → super; then _parent.report_param(name, value)   // store on sub's progress AND mirror to parent's progress, so both are queryable
    report_factor(input, f)           → _parent.report_factor(input, f)         // factors are global; only the parent owns the factors.txt file
    progress_save()                   → super; then _parent.progress_save()     // sub's _file_progress is usually null (no-op); parent does the real .param write
    heartbeat()                       → _parent.heartbeat()                     // BOINC heartbeat is process-wide, only the outer Logging is wired to it
};
```

Things to keep straight:

- **The constructor wires the parent, but the generic production users detach it.** The `SubLogging` ctor does call `set_parent(&parent.progress())` (`logging.h:111`), so a SubLogging *left attached* (e.g. the test harness at `testing.cpp:269`) walks up to its parent. But both generic production users immediately sever the link: `MorrisonGeneric` and `PocklingtonGeneric` call `_logging->progress().set_parent(nullptr)` right after constructing the SubLogging (`morrison.cpp:456`, `pocklington.cpp:250`), and then re-null per inner node via `_logging->progress() = Progress();` (`morrison.cpp:684`, `pocklington.cpp:464`; a default `Progress` has `_parent == nullptr`, `logging.h:46`). So the observation "parent set to nullptr" is a **correct** factual reading of those call sites — but it does not explain lost time: the cause is the wiped `_time_stage`, not the null parent. The fix lives elsewhere; see §9. (Sibling docs `run-hierarchy.md` and `lay-of-the-land.md` describe it this way.)
- **A SubLogging has its own Progress** with its own `_costs`, `_time_total`, etc. *If still attached*, its `update()` updates its own fields *and* walks up to the parent. In `MorrisonGeneric`/`PocklingtonGeneric` the parent is detached, so the sub's `update()` touches only the sub's own fields. Either way, reading a SubLogging's `time_total()` gives the sub's own elapsed; reading the parent's `time_total()` gives whatever the parent has credited durably plus what descendants (if any are attached) have added in the current `_time_stage` window.
- **`progress_save()` cascades.** When the inner code does `_logging->progress_save()`, both the inner SubLogging's `_file_progress` (often null — sub doesn't usually own one) and the parent's persistent `.param` file get written.
- **Default sub level is `LEVEL_WARNING`.** Inner-task chatter is suppressed unless the caller bumps the level. Result lines from the sub still bubble up through `report()`.
- **Prefix applies to forwarded messages.** `_parent.report(prefix+message, level)`. The sub adds its prefix, the parent prints it. Don't set the same prefix on both or you'll get it doubled.

## 7. The driver: how `Task` exercises Progress

`Task::on_state()` (`framework/task.cpp:163`) is the periodic callback inside any long-running task. It calls `_logging->progress().update(progress(), ops())` under three conditions:

1. **Disk-save tick** (every `DISK_WRITE_TIME` seconds, default 300, or on abort):
   `update(progress, ops)` → `state_save()` → `write_state()` → `update(progress, ops)` again → `progress_save()`.
2. **Progress tick** (every `PROGRESS_TIME` seconds, default 10):
   `update(progress, ops)` → `report_progress()`.
3. **Always at the end:** `_logging->heartbeat()` (BOINC uses this).

`InputTask::init()` calls `update(0, ops())` at the start; `InputTask::done()` calls `update(1, ops())` at the end. Together these bracket the task's lifetime so — **when the Progress still has a parent attached** — the parent walk gets a clean elapsed=0 baseline at start and a final settle at the end.

Net: every long task automatically drives `update()` on its own Logging's Progress. If that Logging is a SubLogging **whose parent is still attached** (e.g. the `-test` harness, `testing.cpp:269`), every drive-by also walks up to the outer's `_time_stage`. **But the two generic production users detach it:** `MorrisonGeneric`/`PocklingtonGeneric` call `set_parent(nullptr)` right after constructing the SubLogging (`morrison.cpp:456`, `pocklington.cpp:250`) and re-null the sub's Progress per inner node, so there the drive-by touches only the sub's own fields and never walks up — the outer's durable time is credited by an explicit `update()` on the *outer* Progress instead (see §6 and §9). **For the attached case you don't need to push time anywhere yourself** — the existing wiring already does it. What you do need is for the outer's own `update()` to fire periodically too, so `_time_total` accumulates durably and isn't all riding on the doomed `_time_stage`.

## 8. The persistence model

The only thing persisted across runs is the `.param` text file. Format: lines of `key=value`. The key `time_total` is special — `file_progress()` feeds it to `time_init()` so wall-clock continues from where it left off.

```
time_total=143.872
P=5
next_fft=1
...
```

`progress_save()` writes empty-skipping. So a value of `""` is removed on next save; a value of `"0"` is kept. To clear a param, set it to `""`.

There is no on-disk schema for `_costs` or `_cur_stage`. **Stages are reconstructed in code on each run.** A test class's constructor (or `run()` setup phase) re-builds the same `add_stage(...)` calls deterministically so progress percentages match across resumes. If you rearrange or add stages, percentages from the old `.param` will look wrong until the run completes.

`time_init()` re-anchors `_timer` to "now" and sets `_time_total` to the saved elapsed. `_time_stage` is implicitly zero on a fresh load.

## 9. Common patterns in this codebase

### Outer test class with no sub-tasks (e.g. `Fermat::run`)

```
ctor:                progress.add_stage(task->cost());    // register expected work
run() body:          task->run();                         // Task::on_state drives update() on logging.progress()
result print:        time_total() returns the real elapsed time durably accumulated in _time_total.
```

No extra Progress glue needed. Time accounting just works.

### Outer test that hands inner tasks a SubLogging (e.g. `MorrisonGeneric`, `PocklingtonGeneric`)

```
ctor:
    _logging = make SubLogging(parent_logging)            ← ctor wires the parent walk...
    _logging->progress().set_parent(nullptr)              ← ...but it is detached right away
    register outer-level stages on parent_logging.progress() (one per inner factor, etc.)

run() body, per iteration:
    _logging->progress() = Progress()              ← fresh sub Progress (re-nulls _parent)
    _logging->progress().add_stage(cur_task->cost())
    cur_task->init(..., _logging.get(), ...)       ← inner uses the SubLogging
    cur_task->run()                                ← drives _logging->progress().update()
                                                     on the DETACHED sub; touches only the
                                                     sub's own fields, does NOT walk up
    _logging->progress().next_stage()
    logging.progress().update(progress_stage, 0)   ← THIS LINE, a direct call on the OUTER
                                                     Progress, captures the inner-task elapsed
                                                     time durably into outer._time_total and
                                                     resets outer._time_stage to 0.
    cur_task.reset()
```

That trailing `update()` on the *outer* Progress (`morrison.cpp:732`, `pocklington.cpp:503`) is the upstream fix from `d14b2f7`. Because the sub's Progress is detached, the inner task never adds elapsed time to the outer's `_time_stage`; the outer's durable time is credited solely by this explicit direct call right after the inner task. (Note the *other* outer-side `update()` — the per-iteration heartbeat at `morrison.cpp:490` — only wipes `_time_stage` to 0; if the d14b2f7 line were absent, the inner-task time would simply never be credited at all.)

### Resuming from a checkpoint

`prst.cpp:332` opens the `.param` file, `logging.file_progress(&file_progress)` loads it. `time_init()` sets `_timer` and `_time_total`. From there everything proceeds as a fresh run, except the result line will report the cumulative elapsed time across all sessions.

## 10. Pitfalls — and how to avoid them

### Pitfall A: assuming `_time_stage` persists past the next `update()`

It doesn't. If you measure, store, or compare `_time_stage` across multiple `update()` calls on the same Progress, you have a bug. Anything that needs to persist must land in `_time_total`.

### Pitfall B: thinking the parent walk credits `_time_total`

It only writes to ancestors' `_time_stage`. Their `_time_total` is updated only when *they themselves* call `update()`. If an outer Progress never gets its own `update()`, its `_time_total` stays at whatever `time_init()` set.

### Pitfall C: extending the API instead of using `update()`

A tempting mistake is to add `Progress::add_time(double)` so the caller could push `task->timer()` into `_time_total` directly. But the existing `update()` already captures wall-clock into `_time_total` via the `_timer` anchor. The minimal fix is one outer `update()` call after the inner task — no new API.

Heuristic: **if you're proposing a new method on `Progress`, re-read the existing six methods first.** If your proposal is functionally a recombination of `update()` + `set_parent()` + `add_stage()`, you don't need a new method.

### Pitfall D: calling `update()` on a SubLogging when you meant the parent

Easy to do because the syntax is the same. `_logging->progress().update(...)` only walks UP from the sub; it doesn't write to the parent's `_time_total`. If you want the parent to credit elapsed time durably, call `logging.progress().update(...)` (the outer one) — that's an entirely different code path.

### Pitfall E: re-ordering `add_stage` calls between runs

`.param` doesn't store stage layout. The code is the schema. Add new stages at the end of an existing layout, and ideally only when you need them — re-ordering changes the meaning of `progress_total()` for in-flight checkpoints.

### Pitfall F: setting an empty param to clear it

`progress_save()` skips empties. If you want to clear a param, set its value to `""` and `progress_save()`. Setting it to `"0"` keeps it.

### Pitfall G: relying on `time_op()` for short tasks

`_time_op` is only updated when `op_count` strictly increases. Tasks that pass `op_count = 0` to `update()` (which most outer test-level updates do) leave `_time_op` stale. Don't read `time_op()` from outside `Task::on_state()`-driven code paths.

## 11. Quick reference — when to call what

| You want to…                                                  | Call                                                                            |
|---------------------------------------------------------------|---------------------------------------------------------------------------------|
| Register expected work for a new stage                        | `progress().add_stage(cost)`                                                    |
| Move to the next stage and tick the parent                    | `progress().next_stage()`                                                       |
| Move to the next stage WITHOUT pinging parents                | `progress().skip_stage()`                                                       |
| Periodic heartbeat from inside a long task                    | `_logging->progress().update(progress(), ops())` (already done by Task)         |
| Capture inner-task elapsed time durably to the outer's total  | `logging.progress().update(logging.progress().progress_stage(), 0)`             |
| Persist progress to disk                                      | `logging.progress_save()`                                                       |
| Reload progress on startup                                    | `logging.file_progress(&file_progress)`                                         |
| Re-anchor wall-clock after manual time mutation               | `progress().time_init(progress().time_total())`                                 |
| Emit a user-facing result line                                | `logging.result(success, "...")` + `logging.result_save("...")`                 |
| Emit a one-shot info line                                     | `logging.info("...")`                                                           |
| Emit an overwriting progress line                             | `logging.report_progress()` (auto-formatted) or `report(buf, LEVEL_PROGRESS)`   |

## 12. Open questions / things this doc punts on

- **`heartbeat()` semantics.** Base class is a no-op; BOINC override does real work. Worth its own short doc once the BOINC path is in play.
- **`state_save_flag` / `state_save`.** Used by Task::on_state to coordinate disk writes. Driven by BOINC and net.cpp paths; a future doc on those paths should cover it.
- **Thread-safety.** None of `Logging`, `SubLogging`, or `Progress` is documented as MT-safe, and that's fine under the current architecture. PRST does use multiple threads — `-t <threads>` and `-spin <threads>` are user-facing — but the parallelism lives *inside* GWnum, parallelizing a single bignum multiplication. Above the GWnum layer, exactly one `Task` runs at a time on the main thread, and Logging/Progress are only mutated from that thread. The two known exceptions are benign: `Task::abort()` from the SIGTERM handler only flips a `bool _abort_flag`, and BOINC's heartbeat similarly only nudges a flag. If multi-`Task` parallelism is ever added above GWnum (e.g. parallel inner factors in `MorrisonGeneric`), this whole subsystem needs locking — `update()`'s parent walk in particular is racy by construction.
- **Why `_time_stage` exists at all.** It looks like a vestigial design from a moment when ancestor reporting needed to know "time since I was last asked." A simpler model would be "every Progress's `_time_total` is sole truth; ancestors aggregate on read," but changing that would break the existing reporting semantics.
