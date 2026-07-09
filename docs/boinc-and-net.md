# BOINC & NETPRST — deep dive

PRST has two optional "distributed-computing" front-ends that wrap the same test machinery in a different I/O and lifecycle shell:

- **BOINC** (`boinc.{h,cpp}`) — runs PRST as a BOINC science app: progress and checkpoints go through the BOINC client (via a thin `bow/` wrapper — the `bow_*` calls are all BOINC client-integration hooks: init, event-poll, progress, trickle, finish), and the client can suspend/abort/quit the task at any checkpoint boundary.
- **NETPRST** (`net.{h,cpp}`) — runs PRST as a thin HTTP worker: it polls a coordinator server for a task, ships checkpoints to/from the server instead of local disk, and POSTs the result back.

Both are **compile-time flavors**, gated by `#define`s in `version.h` (both **commented out by default**, `version.h:3-4`). When enabled, each adds a runtime subcommand to `main()`'s parser that `exit()`s into a sibling `main`: `-boinc` → `boinc_main` (`prst.cpp:200`), `-net` → `net_main` (`prst.cpp:203`). Neither replaces `main()`; both re-use `Fermat`/`Pocklington`/`Proof` and the `Task` layer underneath.

> **Read this first.** `boinc.cpp` tracks the current `Run`/`Proof` interface (`run->run(gwstate, file_checkpoint, file_recoverypoint, logging)`, `boinc.cpp:340`). **`net.cpp` does not** — it calls `proof->run(input, gwstate, …)` (`net.cpp:574`), `fermat->run(input, gwstate, …)` (`net.cpp:585`), and `fermat->run(input, gwstate, …, nullptr)` (`net.cpp:591`), each prepending an `input` argument that no current `Fermat::run`/`Proof::run` overload accepts (`fermat.h:23-24`, `proof.h:85,94`). NETPRST predates the run() refactor (`4764aa9` base class + `fe25a5a` run refactor, which moved `input` to a `Run` member and dropped it from the signature) and **will not compile against the present headers**; because `NETPRST` is off, no build ever catches it. Treat `net.cpp` as a design reference, not working code, until it's reconciled (§7, §8).

Source files:
- `src/boinc.h`, `src/boinc.cpp` (build flag `BOINC`; depends on `bow/bow.h`)
- `src/net.h`, `src/net.cpp` (build flag `NETPRST`; depends on `restc-cpp` + Boost)
- `src/version.h:3-4` (the flags), `src/prst.cpp:200,203` (the dispatch)

Prereqs / companions: `logging-and-progress.md` (`BoincLogging`/`NetLogging` are `Logging` subclasses overriding `report`/`report_progress`/`progress_save`/`state_save_flag`/`heartbeat`), `task-lifecycle.md` (`state_save_flag` and `heartbeat` are the `Task`→`Logging` callbacks polled on the checkpoint cadence; `DISK_WRITE_TIME`/`PROGRESS_TIME` govern frequency), `run-hierarchy.md` (both entry points construct `Run` subclasses, but **inline a subset of `Run::create`'s dispatch** — see §3), `proof-system.md` (both wire `-proof save/build/cert`).

## 1. The two entry points and their override surface

Neither flavor calls `Run::create`. Each hand-rolls a *reduced* dispatcher inline (only Fermat / Proth / Pocklington — **no Morrison, Order, or `*Generic`**). A `c == -1` candidate that would pick `Morrison` under `main()` falls to `Fermat::AUTO` here.

The shared shape: subclass `Logging` to redirect output/progress/checkpoints, then run the normal test. The override surfaces:

| `Logging` virtual | `BoincLogging` (`boinc.cpp`) | `NetLogging` (`net.cpp`) |
|---|---|---|
| `report` | hide most messages; on `LEVEL_RESULT` emit `"Testing complete."`/`"Done."` + `bow_report_progress(1.0)` (`:37-51`) | also append the message to a `"stderr"` `NetFile` that gets uploaded (`:39-45`) |
| `report_progress` | `bow_report_progress(progress().progress_total())` (`:53-59`) | stash `time_op` into the current `PRSTTask` (`:47-51`) |
| `progress_save` | `Logging::progress_save()` + `bow_app_checkpointed(...)` (`:61-65`) | copy progress/time (+ FFT desc/len) into the `PRSTTask` (`:53-62`) |
| `state_save_flag` | poll BOINC events (§4) (`:129-173`) | *(not overridden — uses base)* |
| `heartbeat` | base + maybe send a trickle (`:122-127`) | *(not overridden)* |

BOINC overrides the **control** hooks (`state_save_flag`/`heartbeat`) because the client owns the run's lifecycle; NETPRST overrides the **data** hooks (`report`/`progress_save`) because its state lives on the server, not the controlling client.

Networked-file data model (NETPRST only):

- **`NetFile : File`** (`net.h:62-88`) — a `File` whose bytes live on the server. `read_buffer` lazily `GET`s the file; `commit_writer` queues an async `PUT`; `get_writer`/`free_buffer`/`clear` coordinate with the in-flight upload via a shared buffer + mutex.
- **`LLR2NetFile : NetFile`** (`net.h:90-101`) — same, but rewrites the on-disk buffer to/from the LLR2 checkpoint format (magic/version/type/checksum munging) for compatibility with LLR2 clients.
- **`NetContext`** (`net.h:129-177`) — owns the REST clients (one for GETs, one octet-stream `_putter` for PUTs), the current `PRSTTask`, the upload queue + worker future, and the `NetLogging`.
- **`PRSTTask`** (`net.h:24-58`) — the JSON task descriptor (`id`, `type`, `sk`, `sb`, `n`, `c`, `cyclotomic`, `hex`, `proof`, `count`, `time`, `options` map), wired to JSON by `BOOST_FUSION_ADAPT_STRUCT`.

## 2. `boinc_main` — the central function

`boinc.cpp:175-355`. The skeleton, annotated (the up-to-date flavor):

```cpp
int boinc_main(int argc, char *argv[])
{
    bow_init();                                                       // hand the process to BOINC
    Task::DISK_WRITE_TIME = bow_get_checkpoint_seconds(Task::DISK_WRITE_TIME);  // honor client cadence
    BoincLogging logging;
    options.ProofPointFilename = "prsproof";                          // BOINC-fixed proof filenames
    options.ProofProductFilename = "prsprod";

    Config cnfg; cnfg.ignore("-boinc") /* -t/-fft/-proof/-check/-fermat/-time/-ini */
        .check("-notrickle", logging.send_trickle_messages, false)    // BOINC-only knob
        .default_code([&](const char* p){ if (!input.parse(p)) printf("Unknown option %s.\n", p); })
        .parse_args(argc, argv);

    // filenames namespaced by fingerprint (+ ".<a>" or ".cert"); bow_resolve_filename maps BOINC slots
    // --- inline, REDUCED dispatcher (no Run::create) ---
    if (proof_op == Proof::CERT)                 run.reset(proof.release());
    else if (input.c() == 1 && (input.b() != 2 || log2(input.gk()) >= input.n()) && !options.ForceFermat)
    {
        if (input.is_half_factored()) run.reset(new Pocklington(input, options, logging, proof.get()));
        else                          logging.warning("Not enough factors for Pocklington test. …");
    }
    if (!run) run.reset(new Fermat(options.ForceFermat ? Fermat::FERMAT : Fermat::AUTO, input, options, logging, proof.get()));
    if (proof) { proof->fermat().reset(dynamic_cast<Fermat*>(run.release())); /* wire proof files */ run.reset(proof.release()); }

    input.setup(gwstate);
    try {
        logging.progress().time_init(bow_get_starting_elapsed_time());  // resume elapsed time from BOINC
        run->run(gwstate, file_checkpoint, file_recoverypoint, logging); // CURRENT 4-arg signature
        file_progress.clear();
    } catch (const TaskAbortException&) { failed = true; }
    gwstate.done();

    if (!failed)              bow_finish(PRST_EXIT_NORMAL);    // completed (or exit(0) standalone)
    if (!Task::abort_flag())  bow_finish(PRST_EXIT_FAILURE);   // genuine failure → tell BOINC to abort the job
    return PRST_EXIT_NORMAL;                                   // else: BOINC temporary exit (suspend/quit) — just return
}
```

The three-way exit at the end is the BOINC contract: a clean finish reports the exit code via `bow_finish`; a real failure (not a Ctrl-C/quit) reports `FAILURE` so the server reissues the job; a quit/suspend request just returns so BOINC can resume later from the checkpoint. `bow_finish` does not return in the non-standalone case (it `exit()`s through BOINC).

## 3. Field & method reference

**BOINC — the BOW wrapper calls** (`boinc.cpp`, declared in `bow/bow.h`):

| Call | Role |
|---|---|
| `bow_init()` | attach to the BOINC client (or enter standalone mode). |
| `bow_standalone()` | true when run outside a BOINC client — `BoincLogging` then falls back to normal `Logging` output (`:39,55`). |
| `bow_poll_events()` | bitmask of pending client events: `QUIT_NORMAL` / `QUIT_HBT` (lost client) / `ABORT` / `SUSPENDED`. |
| `bow_report_progress(frac)` | push fraction-done to the client. |
| `bow_app_checkpointed(frac)` | tell the client a checkpoint was just written (paired with `progress_save`). |
| `bow_get_checkpoint_seconds(def)` | the client's preferred checkpoint interval → `Task::DISK_WRITE_TIME`. |
| `bow_get_starting_elapsed_time()` | elapsed time already credited (for resume) → `progress().time_init`. |
| `bow_send_trickle_up(variety, frac)` | send a "still alive / N% done" trickle message. |
| `bow_resolve_filename(in, out)` | map a logical name to the BOINC slot path (used for `-proof pack`, `-proof cert`, `-ini`). |
| `bow_finish(code)` | report completion to the client and exit. |

**BOINC trickles** (`boinc.cpp:67-127`): `TRICKLE_PERIOD = 24h`, `TRICKLE_FIRST_REPORT_DELAY = 10min`. The last-trickle timestamp persists in `trickle_ts.txt` so trickles survive restarts; on first run it's back-dated so the first trickle fires ~10 min in (proving the task started). `heartbeat` checks the clock and sends one when due. `-notrickle` clears `BoincLogging::send_trickle_messages`.

**BOINC `state_save_flag`** (`boinc.cpp:129-173`) — polled by `Task` at each checkpoint boundary; the return/throw drives shutdown:
- `QUIT_NORMAL` / `QUIT_HBT` → `info(...)`, `Task::abort()`, return `true` (checkpoint now, exit cleanly to resume later).
- `ABORT` → `throw TaskAbortException()` (the client wants the job gone).
- `SUSPENDED` → `sleep(1)` loop re-polling until unsuspended (or a more urgent event arrives, via `goto check_again`); logs "Suspending/Resuming" at most ~5 times.

**NETPRST — `NetContext`** (`net.h:129-177`): two `restc_cpp::RestClient`s (`_client` for GETs, `_putter` with `Content-Type: application/octet-stream` for PUTs), the `_task`, an `_upload_queue` + single `_uploadF` worker future + `_upload_mutex` + shared `_upload_buffer`. Key methods: `upload`/`upload_queued`/`upload_cancel`/`upload_wait`/`done`, and `uptime()` (seconds since construction, sent on every request).

**NETPRST — `NetFile` overrides** (`net.cpp:72-222`): `read_buffer` GETs `…/prst/<taskid>/<filename>` (with MD5 verify when `hash`), `commit_writer` stores the buffer + queues an upload, `get_writer`/`free_buffer`/`clear` all take `upload_mutex` and reconcile with any in-flight upload by swapping buffers. `add_child` creates child `NetFile`s with dotted names.

## 4. Lifecycles

**BOINC** (one task per process): `bow_init` → parse → build `Run` (+ optional `Proof` wrap) → `input.setup` → `time_init(bow_get_starting_elapsed_time())` → `run->run(...)` (the `Task` loop polls `state_save_flag`/`heartbeat` on the `DISK_WRITE_TIME`/`PROGRESS_TIME` cadence) → `bow_finish` or return. Checkpoints are ordinary local `File`s (`prst_<fp>.ckpt`/`.rcpt`); BOINC just controls *when* they're written and whether to keep running.

**NETPRST** (`net_main`, `net.cpp:379-644`) — an unbounded work loop:
1. `POST prst/new` (workerID, uptime, version) → deserialize a `PRSTTask` from JSON. On failure: sleep 1 min, retry.
2. Build the `InputNum` from the task: `Hex(…)` / `Phi(c,…)` / `k*n!±c` / `k*n#±c` expressions, or `init(sk, sb, n, c)` for plain `k·b^n+c`, or `read` from a `"number"` `NetFile` when `n == 0` (`:481-497`).
3. Configure `Options` from the task `options` map (`write_time`, `FFT_Increment`, `FFT_Safety`, `support=LLR2`, `strong`, `a`) and `proof_op` from the task `proof` string.
4. Build `Pocklington` or `Fermat::AUTO` (+ optional `Proof`), open checkpoint/recovery `NetFile`s (or `LLR2NetFile`s when `support=LLR2`), and run.
5. `net.upload_wait()` for pending checkpoint PUTs; if the task was aborted (server said timed-out/not-found) `Task::abort_reset()` and loop; else `POST prst/res/<taskid>` with the result (`prime` / `prp` / `prp/<res64>` / cert res64), `time`, `version` — retried up to 10×.

**The async upload path** (`NetContext::upload`, `net.cpp:296-371`): `commit_writer` pushes the file onto `_upload_queue`; a single coroutine drains it, `PUT`ting `…/prst/<taskid>/<filename>` with arguments `md5`, `workerID`, `version`, `uptime`, `progress`, `time`, `time_op` plus all progress `params[...]` (via `NetLogging::Params`) and the file bytes as the body. `HttpAuthenticationException` ("task timed out") and `HttpForbiddenException` ("task not found") set `_task->aborted` and `Task::abort()`; other errors sleep 15 s and retry. So a dropped checkpoint upload silently looks like a slow worker (see Pitfalls).

## 5. The networked checkpoint protocol

The clever part of NETPRST is that `File` is an interface, so the entire test stack writes "checkpoints" without knowing they're HTTP round-trips. The endpoints:

| Operation | HTTP | When |
|---|---|---|
| fetch task | `POST <url>prst/new` | top of each work-loop iteration |
| download checkpoint / number | `GET <url>prst/<taskid>/<filename>` | `NetFile::read_buffer` (lazy, MD5-checked) |
| upload checkpoint | `PUT <url>prst/<taskid>/<filename>` | `NetFile::commit_writer` → async queue |
| post result | `POST <url>prst/res/<taskid>` | after the run completes |

**MD5 integrity** (`net.cpp:163-173`): when `File::hash` is set, a downloaded body is hashed and compared to the `MD5` response header; a mismatch `clear()`s the buffer (treated as "no checkpoint"). Uploads carry the writer's `hash_str()` as the `md5` argument.

**LLR2 compatibility** (`LLR2NetFile`, `net.cpp:224-265`): when the task sets `support=LLR2`, checkpoints are stored in LLR2's on-disk layout. The munging is in the header/trailer bytes — and **byte 4 is the `appid`, not a version** (LLR2 stamps `2`, PRST stamps `4`; the header layout is in the framework's `state-serialization.md`, the munging in `checkpoints.md` §5). On read, an LLR2 file (`MAGIC_NUM`, `byte[4]==2`) is rewritten to PRST's form (`byte[4]=4`, the TYPE byte set, and a `StateValue` iteration at offset 12 decremented); on write, the inverse, plus an LLR2 trailer (a zero word, a 32-bit additive checksum summed from offset 8, and 20 zero words). This lets a NETPRST worker resume a checkpoint produced by an LLR2 client and vice-versa — identical logic to the local-disk `LLR2File` in `support.cpp`.

## 6. Pitfalls

- **`net.cpp` is stale and won't compile.** It calls `proof->run(input, gwstate, …)` (`:574`), `fermat->run(input, gwstate, …)` (`:585`), and `fermat->run(input, gwstate, …, nullptr)` (`:591`) — an interface that predates the run() refactor (`4764aa9`/`fe25a5a`); current signatures take no `input` (`fermat.h:23-24`, `proof.h:85,94`). Because `NETPRST` is `#define`d off, CI never builds it, so the rot is invisible. Anyone enabling NETPRST must first reconcile these call sites (drop the `input` arg; `Pocklington::run`'s `Proof*` overload and the proof-owns-fermat wiring in `boinc.cpp:312-323` are the model).
- **Neither flavor uses `Run::create`.** Both inline a reduced dispatcher covering only Fermat / Proth / Pocklington. A `c == -1` (Morrison), `Order`, or `*Generic` candidate silently runs as a probabilistic `Fermat::AUTO` instead of its deterministic test. If a distributed campaign needs Morrison, the dispatch in `boinc_main`/`net_main` has to be extended (or pointed at `Run::create`).
- **Build flags are off by default** (`version.h:3-4`). The default Windows/Linux/Mac builds are plain PRST; enabling `BOINC`/`NETPRST` also pulls in heavy deps (`bow/`, or `restc-cpp` + Boost) and their build configs.
- **Silent failure is the BOINC/NET failure mode.** A missed BOINC heartbeat just looks like a slow client (the server may reissue the job); a dropped NETPRST checkpoint upload retries every 15 s and otherwise looks like a slow worker. There's no loud error path — diagnose from server-side job state, not a crash.
- **BOINC trickle state is a loose file.** `trickle_ts.txt` lives in the working directory; deleting it re-arms the "first trickle in 10 min" back-dating. It is *not* part of the checkpoint, so it doesn't migrate with `.ckpt`/`.rcpt`.
- **`bow_finish` doesn't return** under a real BOINC client; the trailing `return PRST_EXIT_NORMAL` only executes on the temporary-exit (suspend/quit) path. Don't add cleanup after the `bow_finish` calls expecting it to run.

## 7. Quick reference

| | BOINC | NETPRST |
|---|---|---|
| Build flag | `#define BOINC` (`version.h:3`) | `#define NETPRST` (`version.h:4`) |
| Entry | `boinc_main` (`-boinc`, `prst.cpp:200`) | `net_main` (`-net`, `prst.cpp:203`) |
| Deps | `bow/bow.h` | `restc-cpp`, Boost |
| Logging subclass | `BoincLogging` (control hooks) | `NetLogging` (data hooks) |
| Checkpoint storage | local `File` (client controls cadence) | `NetFile` over HTTP PUT/GET |
| Task source | one task = the CLI input | `POST prst/new` loop |
| Tests supported | Fermat / Proth / Pocklington (+ Proof) | same |
| Current with `run()` refactor? | **yes** | **no — won't compile** |

| You want to… | Where |
|---|---|
| Change checkpoint cadence under BOINC | it's `bow_get_checkpoint_seconds` → `Task::DISK_WRITE_TIME` (`boinc.cpp:186`) |
| Disable trickles | `-notrickle` (`boinc.cpp:261`) |
| Add a server-controlled option | the `net.task()->options` map reads in `net_main` (`net.cpp:502-529`) |
| Change the HTTP endpoints | the `RequestBuilder` URLs in `net.cpp` (`:130,324,456,621`) |
| Support a new checkpoint wire format | subclass `NetFile` like `LLR2NetFile` (`net.cpp:224-265`) |

## 8. Open questions / non-coverage

- **Reconciling `net.cpp` with the current interface.** The exact set of edits to make NETPRST compile (and whether the proof-owns-fermat ownership transfer from `boinc.cpp:312-323` should be mirrored) is unverified — nobody has built it post-refactor. A focused task, not a doc.
- **The `bow/` library.** `bow_init`/`bow_poll_events`/`bow_finish`/… are treated here as the BOINC-integration contract; their implementation lives in the `bow/` submodule and is out of scope.
- **restc-cpp coroutine model.** The `ProcessWithPromiseT` lambdas, `Context::Sleep`, and the single-worker upload future are reproduced from the code; the restc-cpp threading/coroutine semantics themselves are a library concern, not covered.
- **The coordinator server protocol.** This doc describes the client side of `prst/new` / `prst/<id>/<file>` / `prst/res/<id>`; the server's job scheduling, the JSON task schema beyond `PRSTTask`'s fields, and the result semantics are defined server-side and not in this repo.
- **`bow_send_trickle_up` payload semantics.** The `"llr_progress"` variety string and how the BOINC server interprets the trickle fraction are a server-config concern.
