#include <boinc_api.h>
#include <diagnostics.h>

#include <unistd.h>  // dup2

#define BOINC
#include "bow.h"

// #define VERBOSE

static APP_INIT_DATA aid;          // misc task startup info, not accessible by other means
static double initial_cpu_time;    // CPU time used in previous runs
static double checkpoint_offset;   // when app checkpointed (in this session)

void bow_init(void)
{
    BOINC_OPTIONS options;

    if (boinc_file_exists(INIT_DATA_FILE))
    {
        boinc_init_diagnostics(
            BOINC_DIAG_DUMPCALLSTACKENABLED |
            BOINC_DIAG_HEAPCHECKENABLED |
            BOINC_DIAG_TRACETOSTDERR |
//          BOINC_DIAG_REDIRECTSTDOUT |
            BOINC_DIAG_REDIRECTSTDERR
        );

        /* Redirect STDOUT to same log file as STDERR */
        dup2(2, 1); // (STDERR_FILENO, STDOUT_FILENO)
        setbuf(stdout, NULL);
    }
    else
    {
        // Don't redirect output stderr.txt if no Boinc
        boinc_init_diagnostics(
            BOINC_DIAG_DUMPCALLSTACKENABLED |
            BOINC_DIAG_HEAPCHECKENABLED |
            BOINC_DIAG_TRACETOSTDERR
//          BOINC_DIAG_REDIRECTSTDOUT |
//          BOINC_DIAG_REDIRECTSTDERR
        );
    }

    boinc_options_defaults(options);
    options.send_status_msgs      = false;  // We'll send our status manually, when progress is updated
    options.direct_process_action = false;  // We'll handle pause/resume/quit manually
    boinc_init_options(&options);

    // get a copy of full Boinc startup data
    boinc_get_init_data(aid);

    // get CPU time spent in previous sessions (really, until checkpoint from which we'll continue)
    boinc_wu_cpu_time(initial_cpu_time);
}

int bow_get_checkpoint_seconds(int defval)
{
    int cp;

    if (boinc_is_standalone())
        cp = defval;
    else
    {
        cp = aid.checkpoint_period;
        if (cp < 60) cp = 60;        // sanity check, limit to 1 minute
    }

#ifdef VERBOSE
    printf("%s: %d\n", __FUNCTION__, cp);
#endif
    return cp;
}

//
// on error, return false and keep old cpu_time
//
static bool get_child_cpu_time(double& cpu_time)
{
/*
    double t;

#ifdef _WIN32
    // return -1 on error
    if (boinc_process_cpu_time(pid_handle, t) < 0)
        return false;

#elif defined(__linux__)
    // return zero time on error
    t = linux_cpu_time(pid);
    if (t == 0)
        return false;
#elif defined(__APPLE__)
    // There's no easy way to get another process's CPU time in Mac OS X?
    // Report runtime, it's better then nothing.
    t = boinc_elapsed_time();
#else
#error How to get child CPU time for your OS?
#endif
    cpu_time = t;
*/

    // Unlike wrapper, things are simpler here because we need to get time of own process.
    // Use Boinc API function, but note that it do not return errors.
    cpu_time = boinc_worker_thread_cpu_time();
    return true;
}

void bow_report_progress(double fraction_done)
{
    static double session_time;  // keep old time on error

    get_child_cpu_time(session_time);  // unchanged on error
    boinc_report_app_status(initial_cpu_time + session_time,
                            initial_cpu_time + checkpoint_offset,
                            fraction_done
                           );
#ifdef VERBOSE
    printf("%s: done=%f, initial_cpu_time=%f, cpu_time=%f, checkpoint_offset=%f\n", __FUNCTION__, fraction_done, initial_cpu_time, session_time, checkpoint_offset);
#endif
}

void bow_app_checkpointed(double fraction_done)
{
    get_child_cpu_time(checkpoint_offset);
#ifdef VERBOSE
    printf("%s: App checkpointed at %f\n", __FUNCTION__, checkpoint_offset);
#endif
    bow_report_progress(fraction_done);
}

unsigned bow_poll_events(void)
{
    unsigned mask = 0;
    BOINC_STATUS status;

    boinc_get_status(&status);

    if (status.quit_request)   mask |= BOW_EVENT_QUIT_NORMAL;
    if (status.no_heartbeat)   mask |= BOW_EVENT_QUIT_HBT;
    if (status.abort_request)  mask |= BOW_EVENT_ABORT;
    if (status.suspended)      mask |= BOW_EVENT_SUSPENDED;

    return mask;
}

void bow_send_trickle_up(const char *variety, double ratio_done)
{
    char buf[512];
    static double session_cpu_time;  // keep old time on error

    get_child_cpu_time(session_cpu_time);  // unchanged on error
    snprintf(buf, sizeof(buf),
        "<trickle_up>\n"
        "   <progress>%f</progress>\n"
        "   <cputime>%f</cputime>\n"
        "   <runtime>%f</runtime>\n"
        "</trickle_up>\n",
        ratio_done,
        initial_cpu_time + session_cpu_time,
        aid.starting_elapsed_time + boinc_elapsed_time()
    );
    boinc_send_trickle_up((char *)variety, buf);
}

int  bow_standalone(void)
{
    return boinc_is_standalone();
    // return false;  // debug
}

void bow_finish(int code)
{
    boinc_finish(code);
}
