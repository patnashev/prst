
#include <cmath>
#include <string.h>
#include <signal.h>

#include "gwnum.h"
#include "cpuid.h"
#include "arithmetic.h"
#include "exception.h"
#include "cmdline.h"
#include "inputnum.h"
#include "file.h"
#include "logging.h"
#include "task.h"
#include "params.h"
#include "fermat.h"
#include "proof.h"
#include "pocklington.h"
#include "morrison.h"
#include "testing.h"
#include "support.h"
#include "version.h"

#ifdef BOINC
int boinc_main(int argc, char *argv[]);
#endif
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
    //  7 Morrison test params

    GWState gwstate;
    Params params;
    int proof_op = Proof::NO_OP;
    int proof_count = 0;
    std::string proof_cert;
    bool supportLLR2 = false;
    bool force_fermat = false;
    InputNum input;
    int log_level = Logging::LEVEL_WARNING;

    CmdLine()
        .value_number("-t", 0, gwstate.thread_count, 1, 256)
        .value_number("-t", ' ', gwstate.thread_count, 1, 256)
        .value_number("-spin", ' ', gwstate.spin_threads, 1, 256)
        .value_enum("-cpu", ' ', gwstate.instructions, Enum<std::string>().add("SSE2", "SSE2").add("AVX", "AVX").add("FMA3", "FMA3").add("AVX512F", "AVX512F"))
        .value_number("-fft", '+', gwstate.next_fft_count, 0, 5)
        .group("-fft")
            .value_number("+", 0, gwstate.next_fft_count, 0, 5)
            .value_number("safety", ' ', gwstate.safety_margin, -10.0, 10.0)
            .check("generic", gwstate.force_general_mod, true)
            .check("info", gwstate.information_only, true)
            .end()
        .group("-proof")
            .exclusive()
                .ex_case()
                    .value_number("save", ' ', proof_count, 2, 1048576)
                        .prev_check(proof_op, Proof::SAVE)
                .optional()
                    .list("name", ' ', ' ')
                        .value_string(params.ProofPointFilename)
                        .value_string(params.ProofProductFilename)
                        .end()
                    .end()
                .ex_case()
                    .value_number("build", ' ', proof_count, 2, 1048576)
                        .prev_check(proof_op, Proof::BUILD)
                .optional()
                    .list("name", ' ', ' ')
                        .value_string(params.ProofPointFilename)
                        .value_string(params.ProofProductFilename)
                        .value_string(proof_cert)
                        .end()
                    .value_string("security", ' ', params.ProofSecuritySeed)
                    .value_number("roots", ' ', params.RootOfUnitySecurity, 0, 64)
                        .prev_code([&] { if (params.RootOfUnitySecurity.value() == 0) params.RootOfUnityCheck = false; })
                    .end()
                .ex_case()
                    .value_string("cert", ' ', proof_cert)
                        .prev_check(proof_op, Proof::CERT)
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
                .prev_check(params.CheckStrong, true)
            .end()
        .group("-fermat")
            .value_number("a", ' ', params.FermatBase, 2, INT_MAX)
            .end()
            .prev_check(force_fermat, true)
        .group("-time")
            .value_number("write", ' ', Task::DISK_WRITE_TIME, 1, INT_MAX)
            .value_number("progress", ' ', Task::PROGRESS_TIME, 1, INT_MAX)
            .end()
        .group("-support")
            .check("LLR2", supportLLR2, true)
            .end()
        .value_enum("-log", ' ', log_level, Enum<int>().add("debug", Logging::LEVEL_DEBUG).add("info", Logging::LEVEL_INFO).add("warning", Logging::LEVEL_WARNING).add("error", Logging::LEVEL_ERROR))
        .check("-d", log_level, Logging::LEVEL_INFO)
        .check_code("-test", [&] { exit(testing_main(argc, argv)); })
#ifdef BOINC
        .check_code("-boinc", [&] { exit(boinc_main(argc, argv)); })
#endif
#ifdef NETPRST
        .check_code("-net", [&] { exit(net_main(argc, argv)); })
#endif
        .check_code("-v", [&] {
                print_banner();
                exit(0);
            })
        .value_code("-q", 0, [&](const char* param) {
                if (param[0] != '\"' && !isdigit(param[0]))
                    return false;
                if (!input.parse(param))
                    return false;
                return true;
            })
        .default_code([&](const char* param) {
            if (!input.parse(param))
                printf("Unknown option %s.\n", param);
            })
        .parse(argc, argv);

    if (input.empty())
    {
        printf("Usage: PRST {\"K*B^N+C\" | \"N!+C\" | \"N#+C\" | \"N\"} <options>\n");
        printf("Options: [-log {debug | info | warning | error}]\n");
        printf("\t[-t <threads>] [-spin <threads>]\n");
        printf("\t[-time [write <sec>] [progress <sec>]]\n");
        printf("\t[-fft+1] [-fft [+<inc>] [safety <margin>] [info]] [-cpu {SSE2 | AVX | FMA3 | AVX512F}]\n");
        printf("\t-fermat [a <a>] \n");
        printf("\t-proof {save <count> | build <count> [security <seed>] [roots <depth>] | cert {<name> | default}} [name <proof> <product> [{<cert> | default}]]\n");
        printf("\t-check [{near | always| never}] [strong [count <count>] [L <L>]] \n");
        return 0;
    }

    Logging logging(gwstate.information_only && log_level > Logging::LEVEL_INFO ? Logging::LEVEL_INFO : log_level);

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
    std::unique_ptr<Morrison> morrison;

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
    else if (input.c() == -1 && !force_fermat)
    {
        if (input.is_factorized_half())
            morrison.reset(new Morrison(input, params, logging));
        else
        {
            std::string factors;
            for (auto it = input.b_factors().begin(); it != input.b_factors().end(); it++)
                factors += (!factors.empty() ? " * " : "") + it->first.to_string() + (it->second > 1 ? "^" + std::to_string(it->second) : "");
            logging.warning("Not enough factors of b for Morrison test. Factorized part: %s.\n", factors.data());
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
        else if (morrison)
        {
            fingerprint = File::unique_fingerprint(fingerprint, morrison->factors());
            File file_checkpoint("prst_" + std::to_string(gwstate.fingerprint) + ".c", fingerprint);
            File file_params("prst_" + std::to_string(gwstate.fingerprint) + ".p", fingerprint);
            morrison->run(input, gwstate, file_checkpoint, file_params, logging);
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
            fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat->a()));
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
