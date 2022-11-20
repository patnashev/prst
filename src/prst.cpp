
#include <cmath>
#include <string.h>
#include <signal.h>

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
#include "testing.h"
#include "support.h"
#include "version.h"
#ifdef NETPRST
int net_main(int argc, char *argv[]);
#endif

using namespace arithmetic;

void sigterm_handler(int signo)
{
    Task::abort();
    printf("Terminating...\n");
    signal(signo, sigterm_handler);
}

int main(int argc, char *argv[])
{
    signal(SIGTERM, sigterm_handler);
    signal(SIGINT, sigterm_handler);
    setvbuf(stdout, NULL, _IONBF, 0);
#if defined(_MSC_VER) && !defined(_DEBUG)
    _set_error_mode(_OUT_TO_STDERR);
    _set_abort_behavior(0, _CALL_REPORTFAULT);
#endif

    File::FILE_APPID = 4;
    // -1 test metadata
    //  0 number
    //  1 checkpoint
    //  2 strong check checkpoint
    //  3 proof product
    //  4 certificate
    //  5 strong check placeholder
    //  6 proof state

    int i;
    GWState gwstate;
    Params params;
    uint64_t maxMem = 0;
    int proof_op = Proof::NO_OP;
    int proof_count = 0;
    std::string proof_cert;
    bool supportLLR2 = false;
    bool force_fermat = false;
    InputNum input;
    int log_level = Logging::LEVEL_WARNING;

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

            case 'f':
                if (argv[i][2] && isdigit(argv[i][2]))
                    gwstate.known_factors = argv[i] + 2;
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    gwstate.known_factors = argv[i];
                }
                else
                    break;
                continue;

            case 'd':
                if (!argv[i][2])
                    log_level = Logging::LEVEL_INFO;
                else
                    break;
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
                    else if (i < argc - 2 && strcmp(argv[i + 1], "build") == 0)
                    {
                        i += 2;
                        proof_op = Proof::BUILD;
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
                        if (proof_op == Proof::BUILD && i < argc - 1 && argv[i + 1][0] != '-')
                        {
                            i++;
                            proof_cert = argv[i];
                        }
                    }
                    else if (i < argc - 2 && strcmp(argv[i + 1], "security") == 0)
                    {
                        i += 2;
                        params.ProofSecuritySeed = argv[i];
                    }
                    else if (i < argc - 2 && strcmp(argv[i + 1], "roots") == 0)
                    {
                        i += 2;
                        params.RootOfUnitySecurity = atoi(argv[i]);
                        if (params.RootOfUnitySecurity.value() == 0)
                            params.RootOfUnityCheck = false;
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
            else if (i < argc - 1 && strcmp(argv[i], "-support") == 0)
            {
                i++;
                if (strcmp(argv[i], "LLR2") == 0)
                    supportLLR2 = true;
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
            else if (strcmp(argv[i], "-test") == 0)
                return testing_main(argc, argv);
#ifdef NETPRST
            else if (strcmp(argv[i], "-net") == 0)
                return net_main(argc, argv);
#endif
            else if (i < argc - 1 && strcmp(argv[i], "-log") == 0)
            {
                i++;
                if (strcmp(argv[i], "debug") == 0)
                    log_level = Logging::LEVEL_DEBUG;
                if (strcmp(argv[i], "info") == 0)
                    log_level = Logging::LEVEL_INFO;
                if (strcmp(argv[i], "warning") == 0)
                    log_level = Logging::LEVEL_WARNING;
                if (strcmp(argv[i], "error") == 0)
                    log_level = Logging::LEVEL_ERROR;
            }
            else if (strcmp(argv[i], "-v") == 0)
            {
                printf("PRST version " PRST_VERSION "." VERSION_BUILD ", GWnum library version " GWNUM_VERSION);
#ifdef GMP
                GMPArithmetic* gmp = dynamic_cast<GMPArithmetic*>(&GiantsArithmetic::default_arithmetic());
                if (gmp != nullptr)
                    printf(", GMP library version %s", gmp->version().data());
#endif
                printf("\n");
                return 0;
            }
        }
        else
        {
            if (!input.parse(argv[i]))
                printf("Unknown option %s.\n", argv[i]);
        }
    if (input.empty())
    {
        printf("Usage: PRST {\"K*B^N+C\" | \"N!+C\" | \"N#+C\" | \"N\"} <options>\n");
        printf("Options: [-log {debug | info | warning | error}]\n");
        printf("\t[-t <threads>] [-spin <threads>]\n");
        printf("\t[-time [write <sec>] [progress <sec>]]\n");
        printf("\t[-fft+1] [-fft [+<inc>] [safety <margin>]] [-cpu {SSE2 | AVX | FMA3 | AVX512F}]\n");
        printf("\t-fermat [a <a>] \n");
        printf("\t-proof {save <count> | build <count> [security <seed>] [roots <depth>] | cert {<name> | default}} [name <proof> <product> [{<cert> | default}]]\n");
        printf("\t-check [{near | always| never}] [strong [count <count>] [L <L>]] \n");
        return 0;
    }

    Logging logging(log_level);

    if (input.bitlen() < 32)
    {
        Giant num = input.value();
        GWASSERT(num.size() == 1);
        logging.info("Trial division test of %s.\n", input.display_text().data());
        if (is_prime(num.data()[0]))
        {
            logging.result(true, "%s is prime!\n", input.display_text().data());
            logging.result_save(input.input_text() + " is prime!\n");
        }
        else
        {
            logging.result(false, "%s is not prime.\n", input.display_text().data());
            logging.result_save(input.input_text() + " is not prime.\n");
        }
        return 0;
    }

    std::unique_ptr<File> file_proofpoint;
    std::unique_ptr<File> file_proofproduct;
    std::unique_ptr<File> file_cert;
    auto newFile = [&](std::unique_ptr<File>& file, const std::string& filename, uint32_t fingerprint, char type = BaseExp::State::TYPE)
    {
        if (supportLLR2)
            file.reset(new LLR2File(filename, gwstate.fingerprint, type));
        else
            file.reset(new File(filename, fingerprint));
    };

    uint32_t fingerprint = input.fingerprint();
    gwstate.fingerprint = fingerprint;
    newFile(file_cert, !proof_cert.empty() && proof_cert != "default" ? proof_cert : "prst_" + std::to_string(fingerprint) + ".cert", fingerprint, Proof::Certificate::TYPE);
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

    try
    {
        File file_progress("prst_" + std::to_string(gwstate.fingerprint), fingerprint);
        file_progress.hash = false;
        logging.progress_file(&file_progress);

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
            newFile(file_proofpoint, !params.ProofPointFilename.empty() ? params.ProofPointFilename : "prst_" + std::to_string(gwstate.fingerprint) + ".proof", fingerprint);
            newFile(file_proofproduct, !params.ProofProductFilename.empty() ? params.ProofProductFilename : "prst_" + std::to_string(gwstate.fingerprint) + ".prod", fingerprint, Proof::Product::TYPE);
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

        file_progress.clear();
    }
    catch (const TaskAbortException&)
    {
    }

    gwstate.done();

    return 0;
}
