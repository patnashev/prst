
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
#include "container.h"
#include "logging.h"
#include "task.h"
#include "options.h"
#include "fermat.h"
#include "proof.h"
#include "pocklington.h"
#include "morrison.h"
#include "order.h"
#include "testing.h"
#include "batch.h"
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
    //  6 proof checkpoint
    //  8 serialized checkpoint
    //  9 LucasV checkpoint
    // 10 LucasUV checkpoint
    // 11 LucasUV strong check checkpoint

    GWState gwstate;
    Options options;
    int proof_op = Proof::NO_OP;
    int proof_count = 0;
    std::string proof_cert;
    std::string proof_pack;
    bool proof_keep = false;
    std::optional<bool> proof_write;
    bool supportLLR2 = false;
    bool force_fermat = false;
    std::vector<Giant> factors;
    InputNum input;
    InputNum order_a;
    bool show_info = false;
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
            .check("generic", gwstate.force_mod_type, 1)
            .check("info", gwstate.information_only, true)
            .end()
        .group("-proof")
            .exclusive()
                .ex_case()
                    .value_number("save", ' ', proof_count, 2, 1048576)
                        .on_check(proof_op, Proof::SAVE)
                .optional()
                    .list("name", ' ', ' ')
                        .value_string(options.ProofPointFilename)
                        .value_string(options.ProofProductFilename)
                        .end()
                    .value_string("pack", ' ', proof_pack)
                    .check("keep", proof_keep, true)
                    .value_enum("write", ' ', proof_write, Enum<bool>().add("fast", false).add("small", true))
                    .end()
                .ex_case()
                    .value_number("build", ' ', proof_count, 2, 1048576)
                        .on_check(proof_op, Proof::BUILD)
                .optional()
                    .list("name", ' ', ' ')
                        .value_string(options.ProofPointFilename)
                        .value_string(options.ProofProductFilename)
                        .end()
                    .value_string("cert", ' ', proof_cert)
                    .value_string("pack", ' ', proof_pack)
                    .value_string("security", ' ', options.ProofSecuritySeed)
                    .value_number("roots", ' ', options.RootOfUnitySecurity, 0, 64)
                        .on_code([&] { if (options.RootOfUnitySecurity.value() == 0) options.RootOfUnityCheck = false; })
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
                .ex_case().check_code("near", [&] { options.CheckNear = true; options.Check = false; }).end()
                .ex_case().check_code("always", [&] { options.CheckNear = false; options.Check = true; }).end()
                .ex_case().check_code("never", [&] { options.CheckNear = false; options.Check = false; }).end()
                .end()
            .group("strong")
                .check("disable", options.CheckStrong, false)
                .value_number("count", ' ', options.StrongCount, 1, 1048576)
                .value_number("L", ' ', options.StrongL, 1, INT_MAX)
                .value_number("L2", ' ', options.StrongL2, 1, INT_MAX)
                .end()
                //.on_check(options.CheckStrong, true)
            .end()
        .group("-fermat")
            .value_number("a", ' ', options.FermatBase, 2, INT_MAX)
            .end()
            .on_check(force_fermat, true)
        .value_code("-order", ' ', [&](const char* param) { return order_a.parse(param, false) && order_a.value() > 1; })
        .group("-factors")
            .list("list", ' ', ',', false)
                .value_code([&](const char* param) {
                        Giant tmp;
                        tmp = param;
                        if (tmp == 0)
                            return false;
                        if (tmp != 1)
                            input.add_factor(factors.emplace_back(std::move(tmp)));
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
                                input.add_factor(factors.emplace_back(std::move(tmp)));
                        }
                    }
                    return true;
                })
            .check("all", options.AllFactors, true)
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
        .check_code("-batch", [&] { exit(batch_main(argc, argv)); })
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
        .value_code("-ini", ' ', [&](const char* param) {
                File ini_file(param, 0);
                ini_file.read_buffer();
                if (ini_file.buffer().empty())
                    printf("ini file not found: %s.\n", param);
                else
                    cnfg.parse_ini(ini_file);
                return true;
            })
        .check("-i", show_info, true)
        .check("-info", show_info, true)
        .value_code("-q", 0, [&](const char* param) {
                if (param[0] != '\"' && !isdigit(param[0]))
                    return false;
                InputNum tmp;
                if (!tmp.parse(param))
                    return false;
                if (show_info)
                    input.print_info();
                input = std::move(tmp);
                for (auto& f : factors)
                    input.add_factor(f);
                return true;
            })
        .default_code([&](const char* param) {
                InputNum tmp;
                if (!tmp.parse(param))
                    printf("Unknown option %s.\n", param);
                else
                {
                    if (show_info)
                        input.print_info();
                    input = std::move(tmp);
                    for (auto& f : factors)
                        input.add_factor(f);
                }
            })
        .parse_args(argc, argv);

    if (input.empty())
    {
        printf("Usage: PRST {\"[K*]B^N[/D]+C\" | \"N![A]+C\" | \"[p]N#+C\" | \"X\" | \"Phi(3,[-]<number>)\"} <mode> <options>\n");
        printf("Mode: default is primality testing\n");
        printf("\t-v\n");
        printf("\t-info\n");
        printf("\t-test\n");
        printf("\t-batch\n");
        printf("\t-order {<a> | \"<number>\"}\n");
        printf("Options:\n");
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
        printf("\t-check [{near | always| never}] [strong [disable] [count <count>] [L <L>]]\n");
        printf("\t-proof save <count> [name <proof> <product>] [pack <name>] [keep]\n");
        printf("\t-proof build <count> [security <seed>] [roots <depth>] [name <proof> <product>] [pack <name>] [cert <name>] [keep]\n");
        printf("\t-proof cert {<name> | default}\n");
        return 0;
    }
    if (show_info)
        input.print_info();
    if (show_info && !gwstate.information_only)
        return 0;

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
            return 1;
        }
        else
        {
            logging.result(false, "%s is not prime.\n", input.display_text().data());
            logging.result_save(input.input_text() + " is not prime.\n");
            return 0;
        }
    }

    uint32_t fingerprint = input.fingerprint();
    gwstate.fingerprint = fingerprint;
    std::string filename_suffix;
    if (!order_a.empty() && order_a.value() > 1)
        filename_suffix = "." + std::to_string(order_a.fingerprint());
    else if (options.FermatBase)
        filename_suffix = "." + std::to_string(options.FermatBase.value());
    else if (proof_op == Proof::CERT)
        filename_suffix = ".cert";
    std::string filename_prefix = "prst_" + std::to_string(fingerprint);
    File file_progress(filename_prefix + filename_suffix + ".param", fingerprint);
    file_progress.hash = false;
    logging.file_progress(&file_progress);

    auto newFile = [&](std::unique_ptr<File>& file, const std::string& filename, uint32_t fingerprint, char type = BaseExp::StateValue::TYPE)
    {
        if (supportLLR2)
            file.reset(new LLR2File(filename, gwstate.fingerprint, type));
        else
            file.reset(new File(filename, fingerprint));
    };

    std::unique_ptr<File> file_cert;
    newFile(file_cert, !proof_cert.empty() && proof_cert != "default" ? proof_cert : filename_prefix + ".cert", fingerprint, Proof::Certificate::TYPE);
    std::unique_ptr<Proof> proof;
    if (proof_op != Proof::NO_OP)
        proof.reset(new Proof(proof_op, proof_count, input, options, *file_cert, logging));

    std::unique_ptr<Fermat> fermat;
    std::unique_ptr<PocklingtonGeneric> pocklington;
    std::unique_ptr<Morrison> morrison;
    std::unique_ptr<Order> order;

    if (!order_a.empty() && order_a.value() > 1)
    {
        if (proof_op != Proof::NO_OP)
            logging.warning("Proofs are not implemented in order mode.\n");
        if (input.type() != InputNum::KBNC || input.c() != 1 || !input.cofactor().empty())
        {
            logging.error("Order can be computed only for fully factored K*B^N+1 primes.\n");
            return 1;
        }
        order.reset(new Order(order_a, input, options, logging));
        fingerprint = File::unique_fingerprint(fingerprint, std::to_string(order_a.fingerprint()));
    }
    else if (proof_op == Proof::CERT)
    {
        fingerprint = File::unique_fingerprint(fingerprint, file_cert->filename());
    }
    else if ((input.type() == InputNum::FACTORIAL || input.type() == InputNum::PRIMORIAL || (input.type() == InputNum::KBNC && input.bitlen() > 1000 && input.n() < 10)) && input.c() == 1 && !force_fermat && !proof)
    {
        input.factorize_f_p();
        if (input.is_half_factored())
           pocklington.reset(new PocklingtonGeneric(input, options, logging));
        else
        {
            logging.warning("Not enough factors for Pocklington test.\n");
            fermat.reset(new Fermat(Fermat::AUTO, input, options, logging, proof.get()));
        }
    }
    else if (input.type() == InputNum::KBNC && input.c() == 1 && (input.b() != 2 || log2(input.gk()) >= input.n()) && !force_fermat)
    {
        if (input.is_half_factored())
            fermat.reset(new Pocklington(input, options, logging, proof.get()));
        else
        {
            std::string factors;
            for (auto it = input.factors().begin(); it != input.factors().end(); it++)
                factors += (!factors.empty() ? " * " : "") + it->first.to_string() + (it->second > 1 ? "^" + std::to_string(it->second) :  "");
            logging.warning("Not enough factors for Pocklington test. Factored part: %s.\n", factors.data());
            fermat.reset(new Fermat(Fermat::AUTO, input, options, logging, proof.get()));
        }
    }
    else if ((input.type() == InputNum::FACTORIAL || input.type() == InputNum::PRIMORIAL || (input.type() == InputNum::KBNC && input.bitlen() > 1000 && input.n() < 10)) && input.c() == -1 && !force_fermat && !proof)
    {
        input.factorize_f_p();
        if (input.is_half_factored())
            morrison.reset(new MorrisonGeneric(input, options, logging));
        else
        {
            logging.warning("Not enough factors for Morrison test.\n");
            fermat.reset(new Fermat(Fermat::AUTO, input, options, logging, proof.get()));
        }
    }
    else if (input.type() == InputNum::KBNC && input.c() == -1 && !force_fermat)
    {
        if (input.is_half_factored())
        {
            morrison.reset(new Morrison(input, options, logging));
            if (proof_op != Proof::NO_OP)
                logging.warning("Proofs are not implemented in Morrison test. Use -fermat first.\n");
        }
        else
        {
            std::string factors;
            for (auto it = input.factors().begin(); it != input.factors().end(); it++)
                factors += (!factors.empty() ? " * " : "") + it->first.to_string() + (it->second > 1 ? "^" + std::to_string(it->second) : "");
            logging.warning("Not enough factors for Morrison test. Factored part: %s.\n", factors.data());
            fermat.reset(new Fermat(Fermat::AUTO, input, options, logging, proof.get()));
        }
    }
    else
        fermat.reset(new Fermat(force_fermat ? Fermat::FERMAT : Fermat::AUTO, input, options, logging, proof.get()));

    std::unique_ptr<container::FileContainer> proof_container;
    if (!proof_pack.empty())
    {
        proof_container.reset(new container::FileContainer(proof_pack != "default" ? proof_pack : filename_prefix + ".pack"));
        if (proof_container->error() != container::container_error::OK && (proof_op == Proof::BUILD || proof_container->error() != container::container_error::EMPTY))
            logging.warning("Pack %s is corrupted.\n", proof_pack.data());
    }
    std::unique_ptr<FilePacked> file_proofpacked;
    std::unique_ptr<File> file_proofpoint;
    std::unique_ptr<File> file_proofproduct;
    if (fermat && proof)
    {
        fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat->a()) + "." + std::to_string(proof->points()[proof_count].pos));
        newFile(file_proofpoint, !options.ProofPointFilename.empty() ? options.ProofPointFilename : filename_prefix + ".proof", fingerprint);
        if (proof_container)
            file_proofproduct.reset(new FilePacked(!options.ProofProductFilename.empty() ? options.ProofProductFilename : "prod", fingerprint, *proof_container));
        else
            newFile(file_proofproduct, !options.ProofProductFilename.empty() ? options.ProofProductFilename : filename_prefix + ".prod", fingerprint, Proof::Product::TYPE);
        proof->init_files(file_proofpoint.get(), file_proofproduct.get(), file_cert.get());
        if (proof_container)
        {
            file_proofpacked.reset(new FilePacked(!options.ProofPointFilename.empty() ? options.ProofPointFilename : "proof", fingerprint, *proof_container));
            proof->file_points()[0] = file_proofpacked->add_child(std::to_string(0), file_proofpoint->fingerprint());
            proof->file_points()[proof->count()] = file_proofpacked->add_child(std::to_string(proof->count()), file_proofpoint->fingerprint());
        }
        if (proof_write)
            for (int i = 1; i < proof_count; i++)
                fermat->task()->points()[i].value = proof_write.value();
    }
    else if (fermat)
        fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat->a()));
    File file_checkpoint(filename_prefix + filename_suffix + ".ckpt", fingerprint);
    File file_recoverypoint(filename_prefix + filename_suffix + ".rcpt", fingerprint);

    if (gwstate.next_fft_count < logging.progress().param_int("next_fft"))
        gwstate.next_fft_count = logging.progress().param_int("next_fft");
    gwstate.maxmulbyconst = options.maxmulbyconst;
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());

    bool success = false;
    bool failed = false;
    try
    {
        if (order)
            order->run(order_a, options, input, gwstate, file_checkpoint, file_recoverypoint, logging);
        else if (proof_op == Proof::CERT)
            proof->run(input, gwstate, file_checkpoint, file_recoverypoint, logging);
        else if (pocklington)
        {
            pocklington->run(input, gwstate, file_checkpoint, file_recoverypoint, logging);
            success = pocklington->success();
        }
        else if (morrison)
        {
            morrison->run(input, gwstate, file_checkpoint, file_recoverypoint, logging);
            success = morrison->success();
        }
        else if (proof)
        {
            fermat->run(input, gwstate, file_checkpoint, file_recoverypoint, logging, proof.get());
            if (proof_op != Proof::BUILD)
                success = fermat->success();

            if (proof_container)
                proof_container->close();
            if (!proof_keep)
            {
                if (proof_op == Proof::SAVE)
                    for (int i = 1; i < proof->count(); i++)
                        proof->file_points()[i]->clear();
                if (proof_op == Proof::BUILD && proof_container)
                    remove(proof_container->filename().data());
                else if (proof_op == Proof::BUILD)
                {
                    if (!proof->Li())
                        proof->file_points()[0]->clear();
                    proof->file_points()[proof->count()]->clear();
                    for (auto& file : proof->file_products())
                        file->clear();
                }
            }
        }
        else if (fermat)
        {
            fermat->run(input, gwstate, file_checkpoint, file_recoverypoint, logging, nullptr);
            success = fermat->success();
        }

        file_progress.clear();
    }
    catch (const TaskAbortException&)
    {
        if (!gwstate.information_only)
            failed = true;
    }

    gwstate.done();

    return success || failed ? 1 : 0;
}
