
#include <cmath>
#include <string.h>
#include <signal.h>

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
    bool proof_keep = false;
    bool supportLLR2 = false;
    bool force_fermat = false;
    InputNum input;
    std::vector<Giant> factors;
    int log_level = Logging::LEVEL_WARNING;
    std::string log_file;

    Config cnfg;
    cnfg.value_number("-t", 0, gwstate.thread_count, 1, 256)
        .value_number("-t", ' ', gwstate.thread_count, 1, 256)
        .value_number("-spin", ' ', gwstate.spin_threads, 0, 256)
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
                        .on_check(proof_op, Proof::SAVE)
                .optional()
                    .list("name", ' ', ' ')
                        .value_string(params.ProofPointFilename)
                        .value_string(params.ProofProductFilename)
                        .end()
                    .check("keep", proof_keep, true)
                    .end()
                .ex_case()
                    .value_number("build", ' ', proof_count, 2, 1048576)
                        .on_check(proof_op, Proof::BUILD)
                .optional()
                    .list("name", ' ', ' ')
                        .value_string(params.ProofPointFilename)
                        .value_string(params.ProofProductFilename)
                        .end()
                    .value_string("cert", ' ', proof_cert)
                    .value_string("security", ' ', params.ProofSecuritySeed)
                    .value_number("roots", ' ', params.RootOfUnitySecurity, 0, 64)
                        .on_code([&] { if (params.RootOfUnitySecurity.value() == 0) params.RootOfUnityCheck = false; })
                    .check("keep", proof_keep, true)
                    .end()
                .ex_case()
                    .value_string("cert", ' ', proof_cert)
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
        .group("-factors")
            .list("list", ' ', ',', false)
                .value_code([&](const char* param) {
                        Giant tmp;
                        tmp = param;
                        if (tmp == 0)
                            return false;
                        if (tmp != 1)
                            factors.push_back(std::move(tmp));
                        return true;
                    })
                .end()
            .value_code("file", ' ', [&](const char* param) {
                    File factor_file(param, 0);
                    factor_file.read_buffer();
                    if (factor_file.buffer().empty())
                        printf("factor file not found: %s.\n", param);
                    else
                    {
                        std::unique_ptr<TextReader> reader(factor_file.get_textreader());
                        Giant tmp;
                        std::string factor;
                        while (reader->read_textline(factor))
                        {
                            tmp = factor;
                            if (tmp != 0 && tmp != 1)
                                factors.push_back(std::move(tmp));
                        }
                    }
                    return true;
                })
            .check("all", params.AllFactors, true)
            .end()
        .group("-time")
            .value_number("write", ' ', Task::DISK_WRITE_TIME, 1, INT_MAX)
            .value_number("progress", ' ', Task::PROGRESS_TIME, 1, INT_MAX)
            .check("coarse", Task::MULS_PER_STATE_UPDATE, Task::MULS_PER_STATE_UPDATE/10)
            .end()
        .group("-support")
            .check("LLR2", supportLLR2, true)
            .end()
        .group("-log")
            .exclusive()
                .ex_case().check("debug", log_level, Logging::LEVEL_DEBUG).end()
                .ex_case().check("info", log_level, Logging::LEVEL_INFO).end()
                .ex_case().check("warning", log_level, Logging::LEVEL_WARNING).end()
                .ex_case().check("error", log_level, Logging::LEVEL_ERROR).end()
                .end()
            .value_string("file", ' ', log_file)
            .end()
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
        .value_code("-ini", ' ', [&](const char* param) {
                File ini_file(param, 0);
                ini_file.read_buffer();
                if (ini_file.buffer().empty())
                    printf("ini file not found: %s.\n", param);
                else
                    cnfg.parse_ini(ini_file);
                return true;
            })
        .default_code([&](const char* param) {
            if (!input.parse(param))
                printf("Unknown option %s.\n", param);
            })
        .parse_args(argc, argv);

    if (input.empty())
    {
        printf("Usage: PRST {\"K*B^N+C\" | \"N!+C\" | \"N#+C\" | \"N\"} <options>\n");
        printf("Options:\n");
        printf("\t-v\n");
        printf("\t-test\n");
        printf("\t-ini <filename>\n");
        printf("\t-log [{debug | info | warning | error}] [file <filename>]\n");
        printf("\t-time [write <sec>] [progress <sec>] [coarse]\n");
        printf("\t-t <threads>\n");
        printf("\t-spin <threads>\n");
        printf("\t-fft+1\n");
        printf("\t-fft [+<inc>] [safety <margin>] [generic] [info]\n");
        printf("\t-cpu {SSE2 | AVX | FMA3 | AVX512F}\n");
        printf("\t-fermat [a <a>]\n");
        printf("\t-factors [list <factor>,...] [file <filename>] [all]\n");
        printf("\t-check [{near | always| never}] [strong [count <count>] [L <L>]]\n");
        printf("\t-proof save <count> [name <proof> <product>] [keep]\n");
        printf("\t-proof build <count> [security <seed>] [roots <depth>] [name <proof> <product>] [cert <name>] [keep]\n");
        printf("\t-proof cert {<name> | default}\n");
        return 0;
    }
    for (auto& f : factors)
        input.add_factor(f);

    Logging logging(gwstate.information_only && log_level > Logging::LEVEL_INFO ? Logging::LEVEL_INFO : log_level);
    if (!log_file.empty())
        logging.file_log(log_file);

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
        {
            morrison.reset(new Morrison(input, params, logging));
            if (proof_op != Proof::NO_OP)
                logging.warning("Proofs are not implemented in Morrison test. Use -fermat first.\n");
        }
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
        fermat.reset(new Fermat(force_fermat ? Fermat::FERMAT : Fermat::AUTO, input, params, logging, proof.get()));


    gwstate.maxmulbyconst = params.maxmulbyconst;
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());

    try
    {
        File file_progress("prst_" + std::to_string(gwstate.fingerprint), fingerprint);
        file_progress.hash = false;
        logging.file_progress(&file_progress);

        if (proof_op == Proof::CERT)
        {
            fingerprint = File::unique_fingerprint(fingerprint, file_cert->filename());
            File file_checkpoint("prst_" + std::to_string(gwstate.fingerprint) + ".cert.c", fingerprint);
            File file_recoverypoint("prst_" + std::to_string(gwstate.fingerprint) + ".cert.r", fingerprint);
            proof->run(input, gwstate, file_checkpoint, file_recoverypoint, logging);
        }
        else if (morrison)
        {
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

            if (!proof_keep)
            {
                if (proof_op == Proof::SAVE)
                    for (int i = 1; i < proof->count(); i++)
                        proof->file_points()[i]->clear();
                if (proof_op == Proof::BUILD)
                {
                    if (!proof->Li())
                        proof->file_points()[0]->clear();
                    proof->file_points()[proof->count()]->clear();
                    for (auto file : proof->file_products())
                        file->clear();
                }
            }
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
