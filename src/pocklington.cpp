
#include <algorithm>
#include <set>
#include <cmath>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "pocklington.h"
#include "integer.h"
#include "exception.h"

using namespace arithmetic;

int genProthBase(Giant& k, uint32_t n);

Pocklington::Pocklington(InputNum& input, Options& options, Logging& logging, Proof* proof) : Fermat(Fermat::POCKLINGTON, input, options, logging, proof)
{
    if (type() != POCKLINGTON || _a < 0)
        return;

    if (logging.progress().param_int("a") != 0)
    {
        _a = logging.progress().param_int("a");
        options.maxmulbyconst = _a;
    }

    _done = 1;
    _tasks.reserve(input.factors().size());
    for (int i = 0; i < input.factors().size(); i++)
        if (logging.progress().param("factor" + std::to_string(i)).empty())
        {
            _tasks.emplace_back(i);
            Giant& factor = input.factors()[i].first;
            Giant tmp = _task_fermat_simple->exp()/factor;
            if (tmp != 1)
            {
                _tasks.back().taskFactor.reset(new CarefulExp(std::move(tmp)));
                _tasks.back().taskCheck.reset(new CarefulExp(factor));
            }
        }
        else
            _done *= power(input.factors()[i].first, input.factors()[i].second);

    if (options.AllFactors)
        _all_factors = options.AllFactors.value();
}

void Pocklington::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof)
{
    if (type() != POCKLINGTON || _a < 0)
    {
        Fermat::run(input, gwstate, file_checkpoint, file_recoverypoint, logging, proof);
        return;
    }

    std::string sa = std::to_string(_a);
    File* checkpoint = file_checkpoint.add_child(sa, File::unique_fingerprint(file_checkpoint.fingerprint(), sa));
    File* recoverypoint = file_recoverypoint.add_child(sa, File::unique_fingerprint(file_recoverypoint.fingerprint(), sa));

    logging.info("Pocklington test of %s, a = %d, complexity = %d.\n", input.display_text().data(), _a, (int)logging.progress().cost_total());
    Fermat::run(input, gwstate, *checkpoint, *recoverypoint, logging, proof);

    Giant tmp;
    while (!_tasks.empty())
    {
        if (!success())
            return;
        std::vector<Giant> G;
        G.reserve(_tasks.size());
        std::vector<int> factors;
        std::string factors_str;

        for (auto it = _tasks.begin(); it != _tasks.end(); )
        {
            if (it->taskFactor)
            {
                it->taskFactor->init_giant(&input, &gwstate, &logging, std::move(*_task->result()));
                it->taskFactor->run();
                *_task->result() = std::move(it->taskFactor->X0());
                tmp = std::move(*it->taskFactor->result());

                it->taskCheck->init_giant(&input, &gwstate, &logging, std::move(tmp));
                it->taskCheck->run();
                if (*it->taskCheck->result() != 1)
                {
                    logging.warning("Arithmetic error, restarting.");
                    continue;
                }
                tmp = std::move(it->taskCheck->X0());
            }
            else
                tmp = *_task->result();
            if (tmp != 1)
            {
                tmp -= 1;
                G.emplace_back(std::move(tmp));
                _done *= power(input.factors()[it->index].first, input.factors()[it->index].second);
                factors.push_back(it->index);
                factors_str += (!factors_str.empty() ? ", " : "") + input.factors()[it->index].first.to_string();
                it = _tasks.erase(it);
            }
            else
                it++;
        }

        if (G.size() > 0)
        {
            if (G.size() > 1)
            {
                Product taskP(G.begin(), G.end());
                taskP.init(&input, &gwstate, &logging);
                taskP.run();
                tmp = std::move(taskP.result());
            }
            else
                tmp = std::move(G[0]);

            logging.set_prefix(input.display_text() + " ");
            if (factors_str.size() > 30)
                factors_str = factors_str.substr(0, 30) + "...";
            logging.info("Checking gcd with factors {%s}.\n", factors_str.data());
            logging.set_prefix("");
            tmp.gcd(*gwstate.N);
            if (tmp != 1)
            {
                _res64 = tmp.to_res64();
                logging.result(_prime, "%s is not prime. Factor RES64: %s.\n", input.display_text().data(), _res64.data());
                logging.result_save(input.input_text() + " is not prime. Factor RES64: " + _res64 + ".\n");
                break;
            }
            for (auto it = factors.begin(); it != factors.end(); it++)
                logging.report_param("factor" + std::to_string(*it), "done");
        }

        if (_tasks.empty() || (!_all_factors && _done.bitlen()*2 + 10 > gwstate.N->bitlen() &&
            (_done.bitlen()*2 > gwstate.N->bitlen() + 10 || square(_done) > *gwstate.N)))
        {
            _prime = true;
            break;
        }
        else
        {
            if (proof != nullptr)
            {
                logging.error("Pocklington test needs to restart, disable proofs to proceed.\n");
                throw TaskAbortException();
            }

            PrimeIterator primes = PrimeIterator::get();
            for (; *primes <= _a; primes++);
            _a = *primes;

            if (!_task->smooth())
            {
                gwstate.done();
                gwstate.maxmulbyconst = _a;
                input.setup(gwstate);
            }

            sa = std::to_string(_a);
            checkpoint = file_checkpoint.add_child(sa, File::unique_fingerprint(file_checkpoint.fingerprint(), sa));
            recoverypoint = file_recoverypoint.add_child(sa, File::unique_fingerprint(file_recoverypoint.fingerprint(), sa));

            logging.report_param("a", sa);
            logging.progress().add_stage(_task->cost());
            logging.progress().update(0, 0);
            logging.progress_save();

            logging.info("Restarting Pocklington test of %s, a = %d.\n", input.display_text().data(), _a);
            Fermat::run(input, gwstate, *checkpoint, *recoverypoint, logging, nullptr);
        }
    }

    if (_prime)
    {
        logging.result(_prime, "%s is prime!\n", input.display_text().data());
        logging.result_save(input.input_text() + " is prime!\n");
    }

    file_checkpoint.clear(true);
    file_recoverypoint.clear(true);
}

PocklingtonGeneric::PocklingtonGeneric(InputNum& input, Options& options, Logging& logging) : _options(options)
{
    if (logging.progress().param_int("a") != 0)
        _a = logging.progress().param_int("a");
    else
        _a = options.FermatBase ? options.FermatBase.value() : 3;
    options.maxmulbyconst = (_a > GWMULBYCONST_MAX ? 1 : _a);

    std::string& st = logging.progress().param("factors");
    for (std::string::const_iterator it = st.begin(), it_p = st.begin(); it != st.end(); )
    {
        while (it != st.end() && *it != ',')
            it++;
        char* str_end;
        if (it != it_p)
            _done_factors.insert((int)std::strtol(&*it_p, &str_end, 10));
        if (it != st.end())
            it++;
        it_p = it;
    }

    Giant tmp;
    Giant exp;
    if (input.cofactor().empty())
        exp = 1;
    else
        exp = input.cofactor();
    _done = 1;
    Giant tmp_exp;
    tmp_exp = 1;
    Giant tmp_done;
    tmp_done = 1;

    for (int i = 0; i < input.factors().size(); i++)
    {
        auto& f = input.factors()[i];
        if (_done_factors.count(i) != 0)
        {
            tmp = f.first;
            tmp.power(f.second);
            tmp_done *= tmp;
            tmp_exp *= tmp;
        }
        else if (f.second > 1)
        {
            tmp = f.first;
            tmp.power(f.second - 1);
            tmp_exp *= tmp;
        }
        if (tmp_done.size() > 8192)
        {
            _done *= tmp_done;
            tmp_done = 1;
        }
        if (tmp_exp.size() > 8192)
        {
            exp *= tmp_exp;
            tmp_exp = 1;
        }
    }
    if (tmp_done != 1)
        _done *= tmp_done;
    if (tmp_exp != 1)
        exp *= tmp_exp;

    _logging.reset(new SubLogging(logging, logging.level() + 1));
    _logging->progress().set_parent(nullptr);
    logging.progress().add_stage(input.bitlen()*std::log2(input.factors().size()));

    create_tasks(input, logging, exp);
}

void PocklingtonGeneric::create_tasks(InputNum& input, Logging& logging, arithmetic::Giant& exp)
{
    //_logging.reset(new SubLogging(logging, logging.level() + 1));

    std::vector<std::unique_ptr<FactorTree>> factors;
    int i = 0;
    auto it = _done_factors.begin();
    for (; i < input.factors().size() && it != _done_factors.end() && *it == i; it++, i++);
    while (i < input.factors().size())
    {
        factors.emplace_back(new FactorTree(input.factors()[i].first, i));
        for (i++; i < input.factors().size() && it != _done_factors.end() && *it == i; it++, i++);
    }

    _tree.reset(new FactorTree(factors));
    _tree->exp() = std::move(exp);
/*
    bool CheckStrong = _options.CheckStrong ? _options.CheckStrong.value() : true;

    arithmetic::Giant debug_value = input.value() - 1;
    std::vector<Giant> debug_exp;
    debug_exp.emplace_back() = 1;
    double cost = 0;
    std::vector<FactorTree<BaseExp>*> stack;
    stack.push_back(_tree.get());
    while (!stack.empty())
    {
        if (!stack.back()->task() && stack.back()->exp() != 1)
        {
            debug_exp.push_back(debug_exp.back()*stack.back()->exp());
            if (stack.back()->is_factor())
            {
                GWASSERT(debug_exp.back() == debug_value);
            }
            int bitlen = stack.back()->exp().bitlen();
            int checks = _options.StrongCount ? _options.StrongCount.value() : 16;
            checks = bitlen > checks*1000 ? checks : bitlen/1000 + 1;

            if (bitlen < 30)
                stack.back()->task().reset(new CarefulExp(std::move(stack.back()->exp())));
            else if (stack.size() == 1 || (stack.size() == 2 && !stack[0]->task()))
            {
                if (CheckStrong)
                    stack.back()->task().reset(new FastLiCheckExp(std::move(stack.back()->exp()), checks, 0));
                else
                    stack.back()->task().reset(new FastExp(std::move(stack.back()->exp())));
            }
            else
            {
                if (CheckStrong)
                    stack.back()->task().reset(new LiCheckExp(std::move(stack.back()->exp()), checks, 0));
                else
                    stack.back()->task().reset(new SlidingWindowExp(std::move(stack.back()->exp())));
            }
            stack.back()->task()->set_error_check(!_options.CheckNear || _options.CheckNear.value(), _options.Check && _options.Check.value());
            _logging->progress().add_stage(stack.back()->task()->cost());
        }
        if (!stack.back()->is_factor() && !stack.back()->left()->task())
            stack.push_back(stack.back()->left());
        else if (!stack.back()->is_last() && !stack.back()->right()->task())
            stack.push_back(stack.back()->right());
        else
        {
            stack.pop_back();
            debug_exp.pop_back();
        }
    }
    logging.progress().add_stage(_logging->progress().cost_total());
    */
}

void PocklingtonGeneric::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging)
{
    _success = false;
    _res64 = "";
    logging.info("Pocklington test of %s, a = %d.\n", input.display_text().data(), _a);
    if (gwstate.information_only)
        throw TaskAbortException();

    bool CheckStrong = _options.CheckStrong ? _options.CheckStrong.value() : true;
    if (CheckStrong)
        logging.info("%s Gerbicz-Li check enabled.\n", input.display_text().data());

    Giant tmp_done;
    tmp_done = 1;
    double pct_bitlen = gwstate.N->bitlen()/1000.0;
    double last_progress = 0;
    auto progress_factors = [&]() {
        if (logging.level() <= Logging::LEVEL_PROGRESS)
        {
            if (tmp_done != 1)
                _done *= tmp_done;
            tmp_done = 1;
            char buf[50];
            snprintf(buf, 50, "%.1f%% of factors tested. ", std::floor(_done.bitlen()/pct_bitlen)/10.0);
            logging.report(buf, Logging::LEVEL_PROGRESS);
        }
        logging.progress().update((_options.AllFactors && _options.AllFactors.value() ? 1 : 2)*_done.bitlen()/pct_bitlen/1000.0, 0);
        logging.heartbeat();
        last_progress = 0;
    };

    while (_res64.empty() && !_prime)
    {
        logging.set_prefix(input.display_text() + " ");
        if (_done > 1)
            logging.info("Restarting with %.1f%% of factors tested.\n", std::floor(_done.bitlen()/pct_bitlen)/10.0);

        Giant tmp, tmp_exp;
        tmp_exp = 1;
        Giant exp;
        exp = 1;
        std::vector<FactorTree*> stack;
        std::vector<Giant> stack_value;
        int task_num = 0;
        std::vector<int> factors;
        std::string factors_str;
        std::vector<Giant> G;
        auto test_G = [&]() {
            if (G.size() > 1)
            {
                Product taskP(G.begin(), G.end());
                taskP.init(&input, &gwstate, _logging.get());
                taskP.run();
                tmp = std::move(taskP.result());
            }
            else
                tmp = std::move(G[0]);

            if (factors_str.size() < 50)
                logging.debug("Checking gcd with factors {%s}.\n", factors_str.data());
            else
                logging.debug("Checking gcd with %d factors.\n", G.size());
            G.clear();
            factors_str.clear();

            tmp.gcd(*gwstate.N);
            if (tmp != 1)
            {
                logging.set_prefix("");
                _res64 = tmp.to_res64();
                logging.result(_prime, "%s is not prime. Factor RES64: %s.\n", input.display_text().data(), _res64.data());
                logging.result_save(input.input_text() + " is not prime. Factor RES64: " + _res64 + ".\n");
                return true;
            }

            for (auto it = factors.begin(); it != factors.end(); it++)
            {
                logging.progress().param("factors") += "," + std::to_string(*it);
                GWASSERT(_done_factors.count(*it) == 0);
                _done_factors.insert(*it);
            }
            factors.clear();

            if (tmp_done != 1)
                _done *= tmp_done;
            tmp_done = 1;
            if (tmp_exp != 1)
                exp *= tmp_exp;
            tmp_exp = 1;

            if (_done_factors.size() == input.factors().size())
                GWASSERT((input.cofactor().empty() ? _done : _done*input.cofactor()) == input.value() - 1);
            if (_done_factors.size() == input.factors().size() || ((!_options.AllFactors || !_options.AllFactors.value()) && _done.bitlen()*2 + 10 > gwstate.N->bitlen() &&
                (_done.bitlen()*2 > gwstate.N->bitlen() + 10 || square(_done) > *gwstate.N)))
            {
                _prime = true;
                return true;
            }
            return false;
        };

        stack.push_back(_tree.get());
        if (_a > GWMULBYCONST_MAX)
        {
            stack_value.emplace_back();
            stack_value.back() = _a;
        }
        while (!stack.empty())
        {
            if (!stack.back()->exp().empty() && stack.back()->exp() != 1)
            {
                if (stack.back()->index() >= 0)
                    GWASSERT(input.factors()[stack.back()->index()].first == stack.back()->exp());
                task_num++;
                int bitlen = stack.back()->exp().bitlen();
                int checks = _options.StrongCount ? _options.StrongCount.value() : 16;
                checks = bitlen > checks*1000 ? checks : bitlen/1000 + 1;

                std::unique_ptr<BaseExp> cur_task;
                if (bitlen < 32)
                    cur_task.reset(new CarefulExp(std::move(stack.back()->exp())));
                else if (stack_value.empty())
                {
                    if (CheckStrong)
                        cur_task.reset(new FastLiCheckExp(std::move(stack.back()->exp()), checks, 0));
                    else
                        cur_task.reset(new FastExp(std::move(stack.back()->exp())));
                }
                else
                {
                    if (CheckStrong)
                        cur_task.reset(new LiCheckExp(std::move(stack.back()->exp()), checks, 0));
                    else
                        cur_task.reset(new SlidingWindowExp(std::move(stack.back()->exp())));
                }
                stack.back()->exp().arithmetic().free(stack.back()->exp());
                cur_task->set_error_check(!_options.CheckNear || _options.CheckNear.value(), _options.Check && _options.Check.value());
                _logging->progress() = Progress();
                _logging->progress().add_stage(cur_task->cost());
                _logging->set_prefix("stage " + std::to_string(task_num) + " ");

                File* checkpoint = nullptr;
                File* recoverypoint = nullptr;
                if (stack.size() == 1 || (stack.size() == 2 && stack[0]->exp() == 1))
                {
                    std::string file_num = std::to_string(_a) + "." + _done.to_res64() + "." + std::to_string(_done.bitlen());
                    checkpoint = file_checkpoint.add_child(file_num, File::unique_fingerprint(file_checkpoint.fingerprint(), file_num));
                    recoverypoint = file_recoverypoint.add_child(file_num, File::unique_fingerprint(file_recoverypoint.fingerprint(), file_num));
                }
                if (!stack_value.empty())
                {
                    if (CarefulExp* task = dynamic_cast<CarefulExp*>(cur_task.get()))
                        task->init_giant(&input, &gwstate, _logging.get(), std::move(stack_value.back()));
                    if (SlidingWindowExp* task = dynamic_cast<SlidingWindowExp*>(cur_task.get()))
                        task->init(&input, &gwstate, checkpoint, _logging.get(), std::move(stack_value.back()));
                    if (LiCheckExp* task = dynamic_cast<LiCheckExp*>(cur_task.get()))
                        task->init(&input, &gwstate, checkpoint, recoverypoint, _logging.get(), std::move(stack_value.back()));
                }
                else
                {
                    if (CarefulExp* task = dynamic_cast<CarefulExp*>(cur_task.get()))
                        task->init_small(&input, &gwstate, _logging.get(), _a);
                    if (FastExp* task = dynamic_cast<FastExp*>(cur_task.get()))
                        task->init(&input, &gwstate, checkpoint, _logging.get(), _a);
                    if (FastLiCheckExp* task = dynamic_cast<FastLiCheckExp*>(cur_task.get()))
                        task->init(&input, &gwstate, checkpoint, recoverypoint, _logging.get(), _a);
                }

                cur_task->run();
                if (!stack_value.empty())
                    stack_value.back() = std::move(cur_task->X0());
                stack_value.push_back(std::move(*cur_task->result()));
                if (stack.size() == 1)
                    exp = std::move(cur_task->exp());
                last_progress += cur_task->timer();
                _logging->progress().next_stage();
                cur_task.reset();

                if (stack.back()->is_factor())
                {
                    if (stack_value.back() != 1 && !_success)
                    {
                        logging.set_prefix("");
                        _res64 = stack_value.back().to_res64();
                        logging.result(_prime, "%s is not prime. RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), logging.progress().time_total());
                        logging.result_save(input.input_text() + " is not prime. RES64: " + _res64 + ", time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
                        break;
                    }
                    if (stack_value.back() != 1)
                    {
                        logging.warning("Arithmetic error, restarting.");
                        break;
                    }
                    _success = true;
                    int index = stack.back()->index();
                    stack.pop_back();
                    stack_value.pop_back();
                    if (stack_value.back() != 1)
                    {
                        _logging->info("Factor #%d added to gcd.\n", index);
                        G.push_back(std::move(stack_value.back()));
                        G.back() -= 1;
                        factors.push_back(index);

                        auto& f = input.factors()[factors.back()];
                        tmp_done *= power(f.first, f.second);
                        if (tmp_done.size() > 8192)
                        {
                            _done *= tmp_done;
                            tmp_done = 1;
                        }
                        tmp_exp *= f.first;
                        if (tmp_exp.size() > 8192)
                        {
                            exp *= tmp_exp;
                            tmp_exp = 1;
                        }
                        if (factors_str.size() < 50)
                            factors_str += (!factors_str.empty() ? ", " : "") + f.first.to_string();

                        if (G.size() == 100)
                        {
                            if (test_G())
                                break;
                            logging.progress_save();
                        }
                    }
                    else
                        _logging->info("Factor #%d can't be tested with a=%d.\n", index, _a);
                    if (last_progress > Task::PROGRESS_TIME)
                        progress_factors();
                }
            }

            if (!stack.back()->is_factor() && !stack.back()->left()->exp().empty())
                stack.push_back(stack.back()->left());
            else if (!stack.back()->is_last() && !stack.back()->right()->exp().empty())
                stack.push_back(stack.back()->right());
            else
            {
                stack.pop_back();
                if (!stack_value.empty())
                    stack_value.pop_back();
            }
        }

        if (!stack.empty())
            break;
        if (!G.empty() && test_G())
            break;

        PrimeIterator primes = PrimeIterator::get();
        for (; *primes <= _a; primes++);
        if (((_options.AllFactors && _options.AllFactors.value()) || input.factors()[0].second > gwstate.N->bitlen()/100) && _done_factors.count(0) == 0)
            while (kronecker(*primes, *gwstate.N) != -1)
                primes++;
        _a = *primes;
        logging.report_param("a", _a);

        gwstate.done();
        gwstate.maxmulbyconst = (_a > GWMULBYCONST_MAX ? 1 : _a);
        input.setup(gwstate);

        create_tasks(input, logging, exp);

        logging.progress_save();
        file_checkpoint.clear(true);
        file_recoverypoint.clear(true);

        logging.set_prefix("");
        logging.info("Restarting Pocklington test of %s, a = %d.\n", input.display_text().data(), _a);
    }

    logging.progress().next_stage();
    logging.set_prefix("");
    if (_prime)
    {
        logging.info("%s %.1f%% of factors tested.\n", input.display_text().data(), std::floor(_done.bitlen()/pct_bitlen)/10.0);
        logging.result(_prime, "%s is prime! Time: %.1f s.\n", input.display_text().data(), logging.progress().time_total());
        logging.result_save(input.input_text() + " is prime! Time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
    }

    file_checkpoint.clear(true);
    file_recoverypoint.clear(true);
}
