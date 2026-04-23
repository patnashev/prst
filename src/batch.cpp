
#include <list>
#include <deque>
#include <tuple>
#include <map>
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

#include "prst.h"
#include "fermat.h"
#include "proof.h"
#include "pocklington.h"
#include "morrison.h"
#include "order.h"
#include "batch.h"
#include "abc_parser.h"

using namespace arithmetic;

int batch_main(int argc, char *argv[])
{
    GWState gwstate;
    Options options;
    bool show_info = false;
    bool trial_division = false;
    int log_level = Logging::LEVEL_WARNING;
    int log_batch_level = Logging::LEVEL_INFO;
    std::string log_file;
    Task::PROGRESS_TIME = 60;

    std::string batch_name;
    bool stop_error = false;
    bool stop_prime = false;
    int stop_composites = 0;
    bool stop_k_prime = false;

    Config cnfg;
    cnfg.ignore("-batch")
        .value_number("-t", 0, gwstate.thread_count, 1, 256)
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
        .group("-check")
            .exclusive()
                .ex_case().check_code("near", [&] { options.CheckNear = true; options.Check = false; }).end()
                .ex_case().check_code("always", [&] { options.CheckNear = false; options.Check = true; }).end()
                .ex_case().check_code("never", [&] { options.CheckNear = false; options.Check = false; }).end()
                .end()
            .group("strong")
                .check("disable", options.CheckStrong, false)
                .value_number("count", ' ', options.StrongCount, 1, 1048576)
                .end()
                //.on_check(options.CheckStrong, true)
            .end()
        .group("-fermat")
            .value_number("a", ' ', options.FermatBase, 2, INT_MAX)
            .end()
            .on_check(options.ForceFermat, true)
        .value_code("-order", ' ', [&](const char* param) { options.OrderA.reset(new InputNum()); return options.OrderA->parse(param, false) && options.OrderA->value() > 1; })
        .group("-factors")
            .check("all", options.AllFactors, true)
            .end()
        .group("-time")
            .value_number("write", ' ', Task::DISK_WRITE_TIME, 1, INT_MAX)
            .value_number("progress", ' ', Task::PROGRESS_TIME, 1, INT_MAX)
            .check("coarse", Task::MULS_PER_STATE_UPDATE, Task::MULS_PER_STATE_UPDATE/10)
            .end()
        .group("-log")
            .exclusive()
                .ex_case().check("debug", log_level, Logging::LEVEL_DEBUG).end()
                .ex_case().check("info", log_level, Logging::LEVEL_INFO).end()
                .ex_case().check("warning", log_level, Logging::LEVEL_WARNING).end()
                .ex_case().check("error", log_level, Logging::LEVEL_ERROR).end()
                .end()
            .value_enum("batch", ' ', log_batch_level, Enum<int>().add("debug", Logging::LEVEL_DEBUG).add("info", Logging::LEVEL_INFO).add("warning", Logging::LEVEL_WARNING).add("error", Logging::LEVEL_ERROR))
            .value_string("file", ' ', log_file)
            .end()
        .check("-d", log_level, Logging::LEVEL_INFO)
        .check("-info", show_info, true)
        .check("-trial", trial_division, true)
        .group("-stop")
            .group("on")
                .check("error", stop_error, true)
                .check("prime", stop_prime, true)
                .value_number("composites", ' ', stop_composites, 1, INT_MAX)
                .check("kprime", stop_k_prime, true)
                .end()
            .end()
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
                batch_name = param;
            })
        .parse_args(argc, argv);

    if (batch_name.empty())
    {
        printf("Usage: PRST -batch {<file> | stdin} <options>\n");
        printf("Options:\n");
        printf("\t-info\n");
        printf("\t-ini <filename>\n");
        printf("\t-log [{debug | info | warning | error}] [batch {debug | info | warning | error}] [file <filename>]\n");
        printf("\t-time [write <sec>] [progress <sec>] [coarse]\n");
        printf("\t-t <threads>\n");
        printf("\t-spin <threads>\n");
        printf("\t-fft+1\n");
        printf("\t-fft [+<inc>] [safety <margin>] [generic] [info]\n");
        printf("\t-cpu {SSE2 | AVX | FMA3 | AVX512F}\n");
        printf("\t-trial\n");
        printf("\t-fermat [a <a>]\n");
        printf("\t-order {<a> | \"K*B^N+C\"}\n");
        printf("\t-factors all\n");
        printf("\t-check [{near | always| never}] [strong [disable] [count <count>]]\n");
        printf("\t-stop [on error] [on prime] [on composites <count>] [on kprime]\n");
        printf("Batch file formats: raw (one expression per line), ABC, ABCD, ABC2\n");
        return PRST_EXIT_NORMAL;
    }

    Logging logging_batch(log_batch_level);
    if (!log_file.empty())
        logging_batch.file_log(log_file);

    // --- Parse batch file via CandidateSource or stdin fallback ---

    std::unique_ptr<CandidateSource> source;

    if (batch_name != "stdin")
    {
        source = parse_batch_file(batch_name, logging_batch);
        if (!source || source->size() == 0)
        {
            logging_batch.warning("Batch %s is empty or could not be parsed.\n", batch_name.data());
            return PRST_EXIT_FAILURE;
        }
    }

    size_t total = source ? source->size() : 0;

    std::string filename_suffix;
    if (options.OrderA && !options.OrderA->empty() && options.OrderA->value() > 1)
        filename_suffix = "." + std::to_string(options.OrderA->fingerprint());
    else if (options.FermatBase)
        filename_suffix = "." + std::to_string(options.FermatBase.value());
    File batch_progress(batch_name + filename_suffix + ".param", 0);
    batch_progress.hash = false;
    logging_batch.file_progress(&batch_progress);
    int cur = logging_batch.progress().param_int("cur");
    if (cur > 0)
        logging_batch.info("Restarting at %d.\n", cur + 1);

    int primes = logging_batch.progress().param_int("primes");
    int composites = logging_batch.progress().param_int("composites");

    // Per-k tracking: k_value -> whether a prime has been found for this k
    std::map<std::string, bool> k_prime_found;

    bool success = false;
    for (; batch_name == "stdin" || cur < (int)total; cur++)
    {
        logging_batch.report_param("cur", cur);
        logging_batch.progress().update(total > 0 ? cur/(double)total : 0, 0);
        logging_batch.progress_save();
        if (success && stop_prime)
        {
            Task::abort();
            break;
        }
        if (stop_composites > 0 && composites >= stop_composites)
        {
            logging_batch.info("Stopping: %d consecutive composites reached.\n", composites);
            Task::abort();
            break;
        }

        // --- Get the next candidate expression ---
        std::string expression;
        std::string k_value;

        if (batch_name == "stdin")
        {
            double time = logging_batch.progress().time_total();
            std::getline(std::cin, expression);
            logging_batch.progress().time_init(time);
            if (expression.empty() || Task::abort_flag())
                break;
        }
        else
        {
            Candidate cand;
            if (!source->get(cur, cand) || Task::abort_flag())
                break;
            expression = std::move(cand.expression);
            k_value = std::move(cand.k_value);
        }

        // Per-k skip: if stop_k_prime is set and we already found a prime for this k
        if (stop_k_prime && !k_value.empty() && k_prime_found.count(k_value) && k_prime_found[k_value])
        {
            logging_batch.info("%d of %d: %s, skipping (prime already found for k=%s).\n",
                cur + 1, (int)total, expression.data(), k_value.data());
            continue;
        }

        InputNum input;
        InputNum::ParseResult res = input.parse(expression);
        if (!res)
        {
            printf("Error parsing %s, pos %d: %s.\n", expression.data(), res.pos + 1, res.message.data());
            if (stop_error)
            {
                Task::abort();
                break;
            }
            continue;
        }
        if (show_info)
        {
            input.print_info();
            if (!gwstate.information_only)
                continue;
        }
        std::string run_name = std::to_string(cur + 1) + " of " + std::to_string(total) + ": " + input.display_text();

        Logging logging(gwstate.information_only && log_level > Logging::LEVEL_INFO ? Logging::LEVEL_INFO : log_level);
        if (!log_file.empty())
            logging.file_log(log_file);
        if (batch_name == "stdin")
            logging.level_result_not_success = Logging::LEVEL_WARNING;

        if (input.bitlen() <= 40)
        {
            if (batch_name != "stdin")
                logging_batch.info("%s, Trial division test.\n", run_name.data());
            else
                logging.info("Trial division test of %s.\n", input.display_text().data());
            auto factors = input.factorize_small();
            GWASSERT(!factors.empty());
            if (factors[0] == input.value() || (factors[0] == 0 && factors.size() == 1))
            {
                logging.result(true, "%s is prime!\n", input.display_text().data());
                logging.result_save(input.input_text() + " is prime!\n");
                continue;
            }
            else
            {
                logging.result(false, "%s is not prime.\n", input.display_text().data());
                logging.result_save(input.input_text() + " is not prime.\n");
                continue;
            }
        }
        else if (trial_division)
        {
            auto factors = input.factorize_small();
            if (!factors.empty())
            {
                if (batch_name != "stdin")
                    logging_batch.info("%s, trial division found factor %d.\n", run_name.data(), factors[0]);
                else
                    logging.info("Trial division of %s found factor %d.\n", input.display_text().data(), factors[0]);
                logging.result(false, "%s is not prime.\n", input.display_text().data());
                logging.result_save(input.input_text() + " is not prime.\n");
                continue;
            }
        }

        options.maxmulbyconst = 1;

        uint32_t fingerprint = input.fingerprint();
        std::string filename_prefix = "prst_" + std::to_string(fingerprint);
        File file_progress(filename_prefix + filename_suffix + ".param", fingerprint);
        file_progress.hash = false;
        logging.file_progress(&file_progress);

        std::unique_ptr<Run> run(Run::create(input, options, logging));
        if (!run)
            continue;
        if (batch_name != "stdin")
            logging_batch.info("%s, %s.\n", run_name.data(), run->name().data());

        fingerprint = run->fingerprint();
        File file_checkpoint(filename_prefix + filename_suffix + ".ckpt", fingerprint);
        File file_recoverypoint(filename_prefix + filename_suffix + ".rcpt", fingerprint);

        GWState gwstate_cur;
        gwstate_cur.copy(gwstate);
        gwstate_cur.information_only = gwstate.information_only;
        if (gwstate_cur.next_fft_count < logging.progress().param_int("next_fft"))
            gwstate_cur.next_fft_count = logging.progress().param_int("next_fft");
        gwstate_cur.maxmulbyconst = options.maxmulbyconst;
        input.setup(gwstate_cur);
        logging.info("Using %s.\n", gwstate_cur.fft_description.data());

        success = false;
        bool failed = false;
        try
        {
            run->run(gwstate_cur, file_checkpoint, file_recoverypoint, logging);
            success = run->success();
            file_progress.clear();
        }
        catch (const TaskAbortException&)
        {
            if (!gwstate.information_only)
                failed = true;
        }

        gwstate_cur.done();

        if (success)
        {
            primes++;
            logging_batch.report_param("primes", primes);
            composites = 0;
            logging_batch.report_param("composites", composites);
            if (!k_value.empty())
                k_prime_found[k_value] = true;
        }
        else if (!failed)
        {
            composites++;
            logging_batch.report_param("composites", composites);
        }
        if (failed && stop_error)
            Task::abort();
        if (Task::abort_flag())
            break;
    }

    logging_batch.progress().update(total > 0 ? cur/(double)total : 0, 0);
    if (Task::abort_flag() && batch_name != "stdin")
    {
        logging_batch.progress_save();
        return PRST_EXIT_FAILURE;
    }
    else
    {
        batch_progress.clear();
        logging_batch.info("Batch of %d, primes: %d, time: %.1f s.\n", cur, primes, logging_batch.progress().time_total());
    }

    return PRST_EXIT_NORMAL;
}
