
#include <cmath>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "arithmetic.h"
#include "exception.h"
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

#define BOINC
#include "bow/bow.h"

using namespace arithmetic;

void BoincLogging::report(const std::string& message, int level)
{
    if (level == LEVEL_RESULT && !bow_standalone())
    {
        bow_report_progress(1.0);  // and hide message
        Logging::report("Testing complete.\n", level);
    }
    else
        Logging::report(message, level);
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
    int i;
    GWState gwstate;
    Params params;
    uint64_t maxMem = 0;
    int proof_op = Proof::NO_OP;
    int proof_count = 0;
    std::string proof_cert;
    bool force_fermat = false;
    InputNum input;

    bow_init();
    Task::DISK_WRITE_TIME = bow_get_checkpoint_seconds(Task::DISK_WRITE_TIME);
    printf("PRST version " PRST_VERSION "." VERSION_BUILD ", GWnum library version " GWNUM_VERSION);  // always print in Boinc mode (data collected by validator)
#ifdef GMP
    GMPArithmetic* gmp = dynamic_cast<GMPArithmetic*>(&GiantsArithmetic::default_arithmetic());
    if (gmp != nullptr)
        printf(", GMP library version %s", gmp->version().data());
#endif
    printf("\n");

    params.ProofPointFilename = "proof";
    params.ProofProductFilename = "prod";

    for (i = 1; i < argc; i++)
        if (argv[i][0] == '-' && argv[i][1])
        {
            switch (argv[i][1])
            {
            case 't':
                if (argv[i][2] && isdigit(argv[i][2]))
                    gwstate.thread_count = atoi(argv[i] + 2);
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    gwstate.thread_count = atoi(argv[i]);
                }
                else
                    break;
                if (gwstate.thread_count == 0 || gwstate.thread_count > 64)
                    gwstate.thread_count = 1;
                continue;

            case 'q':
                if (argv[i][2] != '\"' && !isdigit(argv[i][2]))
                    break;
                if (!input.parse(argv[i] + 2))
                {
                    printf("Invalid number format.\n");
                    return 1;
                }
                continue;
            }

            if (i < argc - 1 && strcmp(argv[i], "-M") == 0)
            {
                i++;
                maxMem = InputNum::parse_numeral(argv[i]);
            }
            else if (strncmp(argv[i], "-fft", 4) == 0 && ((!argv[i][4] && i < argc - 1) || argv[i][4] == '+'))
            {
                if (argv[i][4] == '+')
                    gwstate.next_fft_count = atoi(argv[i] + 5);
                else
                    while (true)
                        if (i < argc - 1 && argv[i + 1][0] == '+')
                        {
                            i++;
                            gwstate.next_fft_count = atoi(argv[i] + 1);
                        }
                        else if (i < argc - 2 && strcmp(argv[i + 1], "safety") == 0)
                        {
                            i += 2;
                            gwstate.safety_margin = atof(argv[i]);
                        }
                        else
                            break;
            }
            else if (strcmp(argv[i], "--nthreads") == 0 && i < argc - 1)  // alias for '-t', set by Boinc
            {
                i++;
                gwstate.thread_count = atoi(argv[i]);
                if (gwstate.thread_count == 0 || gwstate.thread_count > 64)
                    gwstate.thread_count = 1;
            }
            else if (strcmp(argv[i], "-generic") == 0)
                gwstate.force_general_mod = true;
            else if (i < argc - 1 && strcmp(argv[i], "-spin") == 0)
            {
                i++;
                gwstate.spin_threads = atoi(argv[i]);
            }
            else if (i < argc - 1 && strcmp(argv[i], "-cpu") == 0)
            {
                while (true)
                    if (i < argc - 1 && strcmp(argv[i + 1], "SSE2") == 0)
                    {
                        i++;
                        gwstate.instructions = "SSE2";
                    }
                    else if (i < argc - 1 && strcmp(argv[i + 1], "AVX") == 0)
                    {
                        i++;
                        gwstate.instructions = "AVX";
                    }
                    else if (i < argc - 1 && strcmp(argv[i + 1], "FMA3") == 0)
                    {
                        i++;
                        gwstate.instructions = "FMA3";
                    }
                    else if (i < argc - 1 && strcmp(argv[i + 1], "AVX512F") == 0)
                    {
                        i++;
                        gwstate.instructions = "AVX512F";
                    }
                    else
                        break;
            }
            else if (i < argc - 2 && strcmp(argv[i], "-proof") == 0)
            {
                while (true)
                    if (i < argc - 2 && strcmp(argv[i + 1], "save") == 0)
                    {
                        i += 2;
                        proof_op = Proof::SAVE;
                        proof_count = atoi(argv[i]);
                    }
                    else if (i < argc - 2 && strcmp(argv[i + 1], "cert") == 0)
                    {
                        i += 2;
                        proof_op = Proof::CERT;
                        proof_cert = argv[i];
                    }
                    else if (i < argc - 3 && strcmp(argv[i + 1], "name") == 0)
                    {
                        i += 2;
                        params.ProofPointFilename = argv[i];
                        i++;
                        params.ProofProductFilename = argv[i];
                    }
                    else
                        break;
            }
            else if (i < argc - 1 && strcmp(argv[i], "-check") == 0)
            {
                while (true)
                    if (i < argc - 1 && strcmp(argv[i + 1], "near") == 0)
                    {
                        i++;
                        params.CheckNear = true;
                        params.Check = false;
                    }
                    else if (i < argc - 1 && strcmp(argv[i + 1], "always") == 0)
                    {
                        i++;
                        params.CheckNear = false;
                        params.Check = true;
                    }
                    else if (i < argc - 1 && strcmp(argv[i + 1], "never") == 0)
                    {
                        i++;
                        params.CheckNear = false;
                        params.Check = false;
                    }
                    else if (i < argc - 1 && strcmp(argv[i + 1], "strong") == 0)
                    {
                        i++;
                        params.CheckStrong = true;
                        if (i < argc - 2 && strcmp(argv[i + 1], "count") == 0)
                        {
                            i += 2;
                            params.StrongCount = atoi(argv[i]);
                        }
                        if (i < argc - 2 && strcmp(argv[i + 1], "L") == 0)
                        {
                            i += 2;
                            params.StrongL = atoi(argv[i]);
                        }
                        if (i < argc - 2 && strcmp(argv[i + 1], "L2") == 0)
                        {
                            i += 2;
                            params.StrongL2 = atoi(argv[i]);
                        }
                    }
                    else
                        break;
            }
            else if (strcmp(argv[i], "-fermat") == 0)
            {
                force_fermat = true;
                if (i < argc - 2 && strcmp(argv[i + 1], "a") == 0)
                {
                    i += 2;
                    params.FermatBase = atoi(argv[i]);
                }
            }
            else if (strcmp(argv[i], "-time") == 0)
            {
                while (true)
                    if (i < argc - 2 && strcmp(argv[i + 1], "write") == 0)
                    {
                        i += 2;
                        Task::DISK_WRITE_TIME = atoi(argv[i]);
                    }
                    else if (i < argc - 2 && strcmp(argv[i + 1], "progress") == 0)
                    {
                        i += 2;
                        Task::PROGRESS_TIME = atoi(argv[i]);
                    }
                    else
                        break;
            }
        }
        else
        {
            if (!input.parse(argv[i]))
                printf("Unknown option %s.\n", argv[i]);
        }
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
