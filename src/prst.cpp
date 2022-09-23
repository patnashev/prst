#define PRST_VERSION "0.1.0"

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
    //guessCpuType();
    //CPU_FLAGS &= ~CPU_FMA3;

    File::FILE_APPID = 4;

    int i, j;
    GWState gwstate;
    GWArithmetic gw(gwstate);
    Params params;
    uint64_t maxMem = 0;
    int proof_op = Proof::NO_OP;
    int proof_count = 0;
    std::string proof_cert;
    InputNum input;
    std::string toFile;
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
                else if (argv[i + 1][0] == '+')
                {
                    i++;
                    gwstate.next_fft_count = atoi(argv[i] + 1);
                }
            }
            else if (strcmp(argv[i], "-generic") == 0)
                gwstate.force_general_mod = true;
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
                    }
                    else if (i < argc - 2 && strcmp(argv[i + 1], "check") == 0)
                    {
                        i += 2;
                        int check = atoi(argv[i]);
                        if (check > proof_count)
                            params.ProofChecksPerPoint = check/proof_count;
                        else
                            params.ProofPointsPerCheck = proof_count/check;
                    }
                    else if (i < argc - 2 && strcmp(argv[i + 1], "security") == 0)
                    {
                        i += 2;
                        params.ProofSecuritySeed = argv[i];
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
                    else if (i < argc - 1 && strcmp(argv[i + 1], "Gerbicz") == 0)
                    {
                        i++;
                        params.CheckGerbicz = true;
                        if (i < argc - 2 && strcmp(argv[i + 1], "count") == 0)
                        {
                            i += 2;
                            params.GerbiczCount = atoi(argv[i]);
                        }
                        if (i < argc - 2 && strcmp(argv[i + 1], "L") == 0)
                        {
                            i += 2;
                            params.GerbiczL = atoi(argv[i]);
                        }
                        if (i < argc - 2 && strcmp(argv[i + 1], "L2") == 0)
                        {
                            i += 2;
                            params.GerbiczL2 = atoi(argv[i]);
                        }
                    }
                    else
                        break;
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
#ifdef NETPRST
            else if (strcmp(argv[i], "-net") == 0)
                return net_main(argc, argv);
#endif
            else if (i < argc - 1 && strcmp(argv[i], "-file") == 0)
            {
                i++;
                toFile = argv[i];
            }
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
                printf("PRST version " PRST_VERSION ", Gwnum library version " GWNUM_VERSION "\n");
                return 0;
            }
        }
        else
        {
            if (!input.parse(argv[i]))
            {
                File file(argv[i], 0);
                if (!input.read(file))
                {
                    printf("File %s is missing or corrupted.\n", argv[i]);
                    return 1;
                }
            }
        }
    if (input.empty())
    {
        printf("Usage: PRST options {\"K*B^N+C\" | file}\n");
        printf("Options: [-t Threads] [-fft+1] -log {debug | info | warning | error}\n");
        return 0;
    }
    if (!toFile.empty())
    {
        File file(toFile, 0);
        input.write(file);
    }

    Logging logging(log_level);

    uint32_t fingerprint = input.fingerprint();
    File file_cert(!proof_cert.empty() && proof_cert != "default" ? proof_cert : "prst_" + std::to_string(fingerprint) + ".cert", fingerprint);
    std::unique_ptr<Proof> proof;
    if (proof_op != Proof::NO_OP)
    {
        if (proof_op == Proof::CERT && (proof_count = Proof::read_cert_power(file_cert)) == 0)
        {
            logging.error("Invalid certificate file.\n");
            return 0;
        }
        proof.reset(new Proof(proof_op, proof_count, input, params, logging));
    }

    std::unique_ptr<Fermat> fermat;
    
    if (proof_op == Proof::CERT)
    {
    }
    else if (proof)
        fermat.reset(new Fermat(input, params, logging, proof.get()));
    else
        fermat.reset(new Fermat(input, params, logging, nullptr));


    gwstate.maxmulbyconst = params.maxmulbyconst;
    input.setup(gwstate);
    
    try
    {
        File file_progress("prst_" + std::to_string(gwstate.fingerprint), fingerprint);
        file_progress.hash = false;
        logging.progress_file(&file_progress);

        if (proof_op == Proof::CERT)
        {
            fingerprint = File::unique_fingerprint(fingerprint, file_cert.filename());
            File file_checkpoint("prst_" + std::to_string(gwstate.fingerprint) + ".cert.c", fingerprint);
            File file_recoverypoint("prst_" + std::to_string(gwstate.fingerprint) + ".cert.r", fingerprint);
            proof->run(input, gwstate, file_cert, file_checkpoint, file_recoverypoint, logging);
        }
        else if (proof)
        {
            fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat->a()) + "." + std::to_string(proof_count) + "." + std::to_string(proof->points()[proof_count]));
            File file_proofpoint(!params.ProofPointFilename.empty() ? params.ProofPointFilename : "prst_" + std::to_string(gwstate.fingerprint) + ".proof", fingerprint);
            File file_proofproduct(!params.ProofProductFilename.empty() ? params.ProofProductFilename : "prst_" + std::to_string(gwstate.fingerprint) + ".prod", fingerprint);
            proof->init_files(&file_proofpoint, &file_proofproduct, &file_cert);

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
