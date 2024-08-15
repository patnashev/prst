
#include <list>
#include <deque>
#include <tuple>
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
#include "options.h"
#include "fermat.h"
#include "proof.h"
#include "pocklington.h"
#include "morrison.h"
#include "order.h"
#include "batch.h"

using namespace arithmetic;

int batch_main(int argc, char *argv[])
{
    GWState gwstate;
    Options options;
    bool force_fermat = false;
    InputNum order_a;
    bool show_info = false;
    int log_level = Logging::LEVEL_WARNING;
    int log_batch_level = Logging::LEVEL_INFO;
    std::string log_file;
    Task::PROGRESS_TIME = 60;

    std::string batch_name;
    bool stop_error = false;
    bool stop_prime = false;

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
            .on_check(force_fermat, true)
        .value_code("-order", ' ', [&](const char* param) { return order_a.parse(param, false) && order_a.value() > 1; })
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
        .group("-stop")
            .group("on")
                .check("error", stop_error, true)
                .check("prime", stop_prime, true)
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
        printf("\t-fermat [a <a>]\n");
        printf("\t-order {<a> | \"K*B^N+C\"}\n");
        printf("\t-factors all\n");
        printf("\t-check [{near | always| never}] [strong [disable] [count <count>]]\n");
        printf("\t-stop [on error] [on prime]\n");
        return 0;
    }

    Logging logging_batch(log_batch_level);
    if (!log_file.empty())
        logging_batch.file_log(log_file);

    std::vector<std::string> batch;

    if (batch_name != "stdin")
    {
        File batch_file(batch_name, 0);
        batch_file.read_buffer();
        std::unique_ptr<TextReader> reader(batch_file.get_textreader());
        std::string st;
        while (reader->read_textline(st))
            batch.push_back(std::move(st));

        if (batch.empty())
        {
            logging_batch.warning("Batch %s is empty.\n", batch_name.data());
            return 0;
        }
    }
    else
        batch.emplace_back();

    std::string filename_suffix;
    if (!order_a.empty() && order_a.value() > 1)
        filename_suffix = "." + std::to_string(order_a.fingerprint());
    else if (options.FermatBase)
        filename_suffix = "." + std::to_string(options.FermatBase.value());
    File batch_progress(batch_name + filename_suffix + ".param", 0);
    batch_progress.hash = false;
    logging_batch.file_progress(&batch_progress);
    int cur = logging_batch.progress().param_int("cur");
    if (cur > 0)
        logging_batch.info("Restarting at %d.\n", cur + 1);

    int primes = logging_batch.progress().param_int("primes");;
    bool success = false;
    for (; cur < batch.size(); cur++)
    {
        logging_batch.report_param("cur", cur);
        logging_batch.progress().update(cur/(double)batch.size(), 0);
        logging_batch.progress_save();
        if (success && stop_prime)
        {
            Task::abort();
            break;
        }
        if (batch_name == "stdin")
        {
            double time = logging_batch.progress().time_total();
            std::getline(std::cin, batch.back());
            logging_batch.progress().time_init(time);
            if (batch.back().empty() || Task::abort_flag())
                break;
            batch.emplace_back();
        }

        InputNum input;
        if (!input.parse(batch[cur]))
        {
            logging_batch.error("Invalid input: %s\n", batch[cur].data());
            if (stop_error)
            {
                Task::abort();
                break;
            }
            continue;
        }
        else if (show_info)
        {
            input.print_info();
            if (!gwstate.information_only)
                continue;
        }
        else if (batch_name != "stdin")
            logging_batch.info("%d of %d: %s", cur + 1, (int)batch.size(), input.display_text().data());

        Logging logging(gwstate.information_only && log_level > Logging::LEVEL_INFO ? Logging::LEVEL_INFO : log_level);
        if (!log_file.empty())
            logging.file_log(log_file);
        if (batch_name == "stdin")
            logging.level_result_not_success = Logging::LEVEL_WARNING;

        if (input.bitlen() < 32)
        {
            Giant num = input.value();
            GWASSERT(num.size() == 1);
            if (batch_name != "stdin")
                logging_batch.info(", Trial division test.\n");
            else
                logging.info("Trial division test of %s.\n", input.display_text().data());
            if (is_prime(num.data()[0]))
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

        options.maxmulbyconst = 1;

        uint32_t fingerprint = input.fingerprint();
        std::string filename_prefix = "prst_" + std::to_string(fingerprint);
        File file_progress(filename_prefix + filename_suffix + ".param", fingerprint);
        file_progress.hash = false;
        logging.file_progress(&file_progress);

        std::unique_ptr<Fermat> fermat;
        std::unique_ptr<PocklingtonGeneric> pocklington;
        std::unique_ptr<Morrison> morrison;
        std::unique_ptr<Order> order;

        if (!order_a.empty() && order_a.value() > 1)
        {
            if (batch_name != "stdin")
                logging_batch.info(", Order.\n");
            if (input.type() != InputNum::KBNC || input.c() != 1 || !input.cofactor().empty())
            {
                logging.error("Order can be computed only for fully factored K*B^N+1 primes.\n");
                continue;
            }
            order.reset(new Order(order_a, input, options, logging));
            fingerprint = File::unique_fingerprint(fingerprint, std::to_string(order_a.fingerprint()));
        }
        else if ((input.type() == InputNum::FACTORIAL || input.type() == InputNum::PRIMORIAL || (input.type() == InputNum::KBNC && input.bitlen() > 1000 && input.n() < 10)) && input.c() == 1 && !force_fermat)
        {
            input.factorize_f_p();
            if (input.is_half_factored())
            {
                if (batch_name != "stdin")
                    logging_batch.info(", generic Pocklington test.\n");
                pocklington.reset(new PocklingtonGeneric(input, options, logging));
            }
            else
            {
                if (batch_name != "stdin")
                    logging_batch.info(", Fermat test.\n");
                logging.warning("Not enough factors for Pocklington test.\n");
                fermat.reset(new Fermat(Fermat::AUTO, input, options, logging, nullptr));
            }
        }
        else if (input.type() == InputNum::KBNC && input.c() == 1 && (input.b() != 2 || log2(input.gk()) >= input.n()) && !force_fermat)
        {
            if (input.is_half_factored())
            {
                if (batch_name != "stdin")
                    logging_batch.info(", Pocklington test.\n");
                fermat.reset(new Pocklington(input, options, logging, nullptr));
            }
            else
            {
                if (batch_name != "stdin")
                    logging_batch.info(", Fermat test.\n");
                std::string factors;
                for (auto it = input.factors().begin(); it != input.factors().end(); it++)
                    factors += (!factors.empty() ? " * " : "") + it->first.to_string() + (it->second > 1 ? "^" + std::to_string(it->second) : "");
                logging.warning("Not enough factors for Pocklington test. Factored part: %s.\n", factors.data());
                fermat.reset(new Fermat(Fermat::AUTO, input, options, logging, nullptr));
            }
        }
        else if ((input.type() == InputNum::FACTORIAL || input.type() == InputNum::PRIMORIAL || (input.type() == InputNum::KBNC && input.bitlen() > 1000 && input.n() < 10)) && input.c() == -1 && !force_fermat)
        {
            input.factorize_f_p();
            if (input.is_half_factored())
            {
                if (batch_name != "stdin")
                    logging_batch.info(", generic Morrison test.\n");
                morrison.reset(new MorrisonGeneric(input, options, logging));
            }
            else
            {
                if (batch_name != "stdin")
                    logging_batch.info(", Fermat test.\n");
                logging.warning("Not enough factors for Morrison test.\n");
                fermat.reset(new Fermat(Fermat::AUTO, input, options, logging, nullptr));
            }
        }
        else if (input.type() == InputNum::KBNC && input.c() == -1 && !force_fermat)
        {
            if (input.is_half_factored())
            {
                morrison.reset(new Morrison(input, options, logging));
                if (batch_name != "stdin" && morrison->is_LLR())
                    logging_batch.info(", Morrison (LLR) test.\n");
                else if (batch_name != "stdin")
                    logging_batch.info(", Morrison test.\n");
            }
            else
            {
                if (batch_name != "stdin")
                    logging_batch.info(", Fermat test.\n");
                std::string factors;
                for (auto it = input.factors().begin(); it != input.factors().end(); it++)
                    factors += (!factors.empty() ? " * " : "") + it->first.to_string() + (it->second > 1 ? "^" + std::to_string(it->second) : "");
                logging.warning("Not enough factors for Morrison test. Factored part: %s.\n", factors.data());
                fermat.reset(new Fermat(Fermat::AUTO, input, options, logging, nullptr));
            }
        }
        else
        {
            fermat.reset(new Fermat(force_fermat ? Fermat::FERMAT : Fermat::AUTO, input, options, logging, nullptr));
            if (batch_name != "stdin" && fermat->type() == Fermat::PROTH)
                logging_batch.info(", Proth test.\n");
            else if (batch_name != "stdin")
                logging_batch.info(", Fermat test.\n");
        }

        if (fermat)
            fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat->a()));
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
            if (order)
                order->run(order_a, options, input, gwstate_cur, file_checkpoint, file_recoverypoint, logging);
            else if (pocklington)
            {
                pocklington->run(input, gwstate_cur, file_checkpoint, file_recoverypoint, logging);
                success = pocklington->success();
            }
            else if (morrison)
            {
                morrison->run(input, gwstate_cur, file_checkpoint, file_recoverypoint, logging);
                success = morrison->success();
            }
            else if (fermat)
            {
                fermat->run(input, gwstate_cur, file_checkpoint, file_recoverypoint, logging, nullptr);
                success = fermat->success();
            }

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
        }
        if (failed && stop_error)
            Task::abort();
        if (Task::abort_flag())
            break;
    }

    logging_batch.progress().update(cur/(double)batch.size(), 0);
    if (Task::abort_flag() && batch_name != "stdin")
        logging_batch.progress_save();
    else
    {
        batch_progress.clear();
        logging_batch.info("Batch of %d, primes: %d, time: %.1f s.\n", cur, primes, logging_batch.progress().time_total());
    }

    return 0;
}
