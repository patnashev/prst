
#include <cmath>
#include <string.h>
#include <iostream>

#include "gwnum.h"
#include "cpuid.h"
#include "arithmetic.h"
#include "exception.h"
#include "config.h"
#include "inputnum.h"
#include "file.h"
#include "logging.h"
#include "task.h"
#include "params.h"
#include "fermat.h"
#include "proof.h"
#include "pocklington.h"
#include "boinc.h"
#include "version.h"

#ifdef _WIN32
#include <windows.h>
#define sleep(x) Sleep(1000*(x))
#else
#include <unistd.h>  // sleep
#endif

#ifdef BOINC

#include "bow/bow.h"

using namespace arithmetic;

void BoincLogging::report(const std::string& message, int level)
{
    if (bow_standalone())
        Logging::report(message, level);
    else if (level == LEVEL_RESULT)
    {
        bow_report_progress(1.0);  // and hide message
        if (progress().cur_stage() <= 1)
            std::cout << "Testing complete.\n";
        else
            std::cout << "Done.\n";
    }
    else
        std::cout << message;
}

void BoincLogging::report_progress()
{
    if (bow_standalone())
        Logging::report_progress();
    else
        bow_report_progress(progress().progress_total());
}

void BoincLogging::progress_save()
{
    bow_app_checkpointed(progress().progress_total());
}

#define TRICKLE_PERIOD             (24 * 3600)   // Send trickles each 24 hours
#define TRICKLE_FIRST_REPORT_DELAY (10 * 60)     // Send first trickle if task is running more then 10 minutes

static const char trickle_file[] = "trickle_ts.txt";

static void save_trickle_file(time_t ts)
{
    FILE *f = fopen(trickle_file, "wt");
    if (f)
    {
        fprintf(f, "%ld\n", (long)ts);
        fclose(f);
    }
}

static void send_trickle_message(Progress &progr)
{
    static time_t last_trickle_time;
    static const char variety[] = "llr_progress";

    time_t now = time(NULL);

    // On first run, try to load saved timestamp of last trickle
    if (last_trickle_time == 0)
    {
        FILE *f = fopen(trickle_file, "rt");
        if (f)
        {
            long tmp;
            if (fscanf(f, "%ld", &tmp) == 1)
                last_trickle_time = tmp;
            fclose(f);
        }
        // If no trickles were sent yet, schedule it to be sent few minutes after start
        // (to be sure that task started up just fine). Otherwise, if Boinc starts a
        // task too close to deadline and did't finish it in time (wrong completion estimate
        // or paused by user), server will be not aware that task is running and will
        // resend potentially good task.
        if (last_trickle_time == 0)
        {
            last_trickle_time = now - TRICKLE_PERIOD + TRICKLE_FIRST_REPORT_DELAY;
            save_trickle_file(last_trickle_time);
        }
    }

    // Time to send new trickle?
    if (now - last_trickle_time >= TRICKLE_PERIOD /* && ratio_done >= 0 */)  /* ratio_done < 0 not applicable here */
    {
        double ratio_done = progr.progress_total();
        bow_send_trickle_up(variety, ratio_done);
        last_trickle_time = now;
        save_trickle_file(last_trickle_time);
    }
}

void BoincLogging::heartbeat()
{
    Logging::heartbeat();
    send_trickle_message(progress());
}

bool BoincLogging::state_save_flag()
{
    if (Task::abort_flag())  return false;   // Don't do anything if already requested to quit

    unsigned mask = bow_poll_events();

check_again:
    if (mask & BOW_EVENT_QUIT_NORMAL)
    {
        error("Exiting - requested by client\n");
        Task::abort();
        return true;    // notice abort ASAP
    }
    else if (mask & BOW_EVENT_QUIT_HBT)
    {
        error("Exiting - lost connection with Boinc client\n");
        Task::abort();
        return true;    // notice abort ASAP
    }
    else if (mask & BOW_EVENT_ABORT)
    {
        error("Boinc requested us to abort\n");
        throw TaskAbortException();
    }
    else if (mask & BOW_EVENT_SUSPENDED)
    {
        static int log_count;

        if (log_count < 5)  // Don't generate too many same messages if user configured Boinc in "suspend when active" mode
            info("Suspending\n");
        do
        {
            sleep(1);
            mask = bow_poll_events();
            if (mask & ~BOW_EVENT_SUSPENDED) goto check_again;  // more critical event received
        } while (mask & BOW_EVENT_SUSPENDED);
        if (log_count < 5)
        {
            info("Resuming\n");
            log_count++;
        }
    }

    return false;
}

int boinc_main(int argc, char *argv[])
{
    GWState gwstate;
    Params params;
    int proof_op = Proof::NO_OP;
    int proof_count = 0;
    std::string proof_cert;
    bool force_fermat = false;
    InputNum input;

    bow_init();
    Task::DISK_WRITE_TIME = bow_get_checkpoint_seconds(Task::DISK_WRITE_TIME);
    print_banner();

    params.ProofPointFilename = "prsproof";
    params.ProofProductFilename = "prsprod";

    Config()
        .ignore("-boinc")
        .value_number("-t", 0, gwstate.thread_count, 1, 256)
        .value_number("-t", ' ', gwstate.thread_count, 1, 256)
        .value_number("--nthreads", ' ', gwstate.thread_count, 1, 256)  // alias for '-t', set by Boinc
        .value_number("-spin", ' ', gwstate.spin_threads, 1, 256)
        .value_enum("-cpu", ' ', gwstate.instructions, Enum<std::string>().add("SSE2", "SSE2").add("AVX", "AVX").add("FMA3", "FMA3").add("AVX512F", "AVX512F"))
        .value_number("-fft", '+', gwstate.next_fft_count, 0, 5)
        .group("-fft")
            .value_number("+", 0, gwstate.next_fft_count, 0, 5)
            .value_number("safety", ' ', gwstate.safety_margin, -10.0, 10.0)
            .check("generic", gwstate.force_general_mod, true)
            .end()
        .group("-proof")
            .exclusive()
                .ex_case()
                    .value_number("save", ' ', proof_count, 2, 1048576)
                        .on_check(proof_op, Proof::SAVE)
                .optional()
                    .list("name", ' ', ' ')
                        .value_string(params.ProofPointFilename)
                        .value_string(params.ProofProductFilename)
                        .end()
                    .end()
                .ex_case()
                    //.value_string("cert", ' ', proof_cert)
                    .value_code("cert", ' ', [&](const char* cert) { bow_resolve_filename(cert, proof_cert); return true; })
                        .on_check(proof_op, Proof::CERT)
                    .end()
                .end()
            .end()
        .group("-check")
            .exclusive()
                .ex_case().check_code("near", [&] { params.CheckNear = true; params.Check = false; }).end()
                .ex_case().check_code("always", [&] { params.CheckNear = false; params.Check = true; }).end()
                .ex_case().check_code("never", [&] { params.CheckNear = false; params.Check = false; }).end()
                .end()
            .group("strong")
                .value_number("count", ' ', params.StrongCount, 1, 1048576)
                .value_number("L", ' ', params.StrongL, 1, INT_MAX)
                .value_number("L2", ' ', params.StrongL2, 1, INT_MAX)
                .end()
                .on_check(params.CheckStrong, true)
            .end()
        .group("-fermat")
            .value_number("a", ' ', params.FermatBase, 2, INT_MAX)
            .end()
            .on_check(force_fermat, true)
        .group("-time")
            .value_number("write", ' ', Task::DISK_WRITE_TIME, 1, INT_MAX)
            .value_number("progress", ' ', Task::PROGRESS_TIME, 1, INT_MAX)
            .end()
        .default_code([&](const char* param) {
            if (!input.parse(param))
                printf("Unknown option %s.\n", param);
            })
        .parse_args(argc, argv);

    if (input.empty())
    {
        printf("No input.\n");
        return 0;
    }

    BoincLogging logging;

    std::unique_ptr<File> file_proofpoint;
    std::unique_ptr<File> file_proofproduct;
    std::unique_ptr<File> file_cert;

    uint32_t fingerprint = input.fingerprint();
    gwstate.fingerprint = fingerprint;
    file_cert.reset(new File(proof_cert, fingerprint));
    std::unique_ptr<Proof> proof;
    if (proof_op != Proof::NO_OP)
        proof.reset(new Proof(proof_op, proof_count, input, params, *file_cert, logging));

    std::unique_ptr<Fermat> fermat;
    
    if (proof_op == Proof::CERT)
    {
    }
    else if (input.c() == 1 && input.b() != 2 && !force_fermat)
    {
        if (input.is_factorized_half())
            fermat.reset(new Pocklington(input, params, logging, proof.get()));
        else
        {
            std::string factors;
            for (auto it = input.b_factors().begin(); it != input.b_factors().end(); it++)
                factors += (!factors.empty() ? " * " : "") + it->first.to_string() + (it->second > 1 ? "^" + std::to_string(it->second) :  "");
            logging.warning("Not enough factors of b for Pocklington test. Factorized part: %s.\n", factors.data());
            fermat.reset(new Fermat(Fermat::AUTO, input, params, logging, proof.get()));
        }
    }
    else
        fermat.reset(new Fermat(Fermat::AUTO, input, params, logging, proof.get()));


    gwstate.maxmulbyconst = params.maxmulbyconst;
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());

    bool failed = false;
    try
    {
        logging.progress().time_init(bow_get_starting_elapsed_time());

        if (proof_op == Proof::CERT)
        {
            fingerprint = File::unique_fingerprint(fingerprint, file_cert->filename());
            File file_checkpoint("prst_" + std::to_string(gwstate.fingerprint) + ".cert.c", fingerprint);
            File file_recoverypoint("prst_" + std::to_string(gwstate.fingerprint) + ".cert.r", fingerprint);
            proof->run(input, gwstate, file_checkpoint, file_recoverypoint, logging);
        }
        else if (proof)
        {
            fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat->a()) + "." + std::to_string(proof->points()[proof_count]));
            file_proofpoint.reset(new File(params.ProofPointFilename, fingerprint));
            file_proofproduct.reset(new File(params.ProofProductFilename, fingerprint));
            proof->init_files(file_proofpoint.get(), file_proofproduct.get(), file_cert.get());

            File file_checkpoint("prst_" + std::to_string(gwstate.fingerprint) + ".c", fingerprint);
            File file_recoverypoint("prst_" + std::to_string(gwstate.fingerprint) + ".r", fingerprint);
            fermat->run(input, gwstate, file_checkpoint, file_recoverypoint, logging, proof.get());
        }
        else if (fermat)
        {
            File file_checkpoint("prst_" + std::to_string(gwstate.fingerprint) + ".c", fingerprint);
            File file_recoverypoint("prst_" + std::to_string(gwstate.fingerprint) + ".r", fingerprint);
            fermat->run(input, gwstate, file_checkpoint, file_recoverypoint, logging, nullptr);
        }
    }
    catch (const TaskAbortException&)
    {
        failed = true;
    }

    gwstate.done();

    if (!failed)             bow_finish(0);   // Boinc task completed, or exit(0) in standalone mode
    if (!Task::abort_flag()) bow_finish(1);   // Failed and it's NOT a Ctrl-C (or quit request), abort Boinc job
    // otherwise it's Boinc temporary exit, just return from program.

    return 0;
}

#endif
