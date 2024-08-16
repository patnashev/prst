
#include <cmath>
#include <algorithm>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "morrison.h"
#include "integer.h"
#include "pocklington.h"

using namespace arithmetic;

// https://eprint.iacr.org/2023/195
// BLS https://doi.org/10.1090/S0025-5718-1975-0384673-1
// For Q = 1  : V_l=2 (mod Factor), l|N+1 for all (D,Factor)
// For Q = -1 : gcd(U_{(N+1)/q}, N) = gcd(V_{(N+1)/2q}, N) due to BLS proof of Theorem 14.
// For Q = -1 factor 2 is tested for free.

Morrison::Morrison(InputNum& input, Options& options, Logging& logging)
{
    int i;
    Giant tmp;
    std::vector<std::pair<Giant, int>> factors;

    bool CheckStrong = options.CheckStrong ? options.CheckStrong.value() : true;
    if (options.AllFactors)
        _all_factors = options.AllFactors.value();

    LucasVMulFast* taskV = nullptr;
    if (!CheckStrong)
        _task.reset(taskV = new LucasVMulFast());

    Giant exp;
    Giant exp_morrison;
    exp = 1;
    exp_morrison = 1;
    int n = 0;

    _factor_tasks.reserve(input.factors().size());
    for (i = 0; i < input.factors().size(); i++)
    {
        auto& factor = input.factors()[i];
        if (factor.first == 2)
        {
            n = factor.second;
            _negQ = (n > 1);
            if (!_negQ)
                _factor_tasks.emplace_back(i);
        }
        else
        {
            if (factor.second > 1)
            {
                Giant divisor;
                if (CheckStrong)
                {
                    divisor = power(factor.first, factor.second - 1);
                    exp_morrison *= divisor;
                    divisor *= factor.first;
                }
                else
                    divisor = power(factor.first, factor.second);
                exp *= divisor;
            }
            else
                exp *= factor.first;
            _factor_tasks.emplace_back(i);
        }
    }
    if (n == 0)
    {
        _task.reset();
        logging.result(false, "%s is not prime, divisible by 2.\n", input.display_text().data());
        logging.result_save(input.input_text() + " is not prime, divisible by 2.\n");
        return;
    }
    if (!input.cofactor().empty())
    {
        exp *= input.cofactor();
        if (CheckStrong)
            exp_morrison *= input.cofactor();
        else
            taskV->mul_giant(input.cofactor(), 1);
    }

    if (log2(exp) < n)
    {
        _LLR = true;
        _factor_tasks.clear();
        if (CheckStrong)
            exp_morrison = exp;
    }
    exp <<= n;
    exp -= 1;
    n--;

    if (CheckStrong)
    {
        exp_morrison <<= n;
        int checks = options.StrongCount ? options.StrongCount.value() : 16;
        _task.reset(new LucasUVMulFast(std::move(exp_morrison), checks));
        options.maxmulbyconst = 2;
    }
    else
        options.maxmulbyconst = 1;

    if (_factor_tasks.size() > 0)
        _taskCheck.reset(new LucasVMulFast(true));
    if (_factor_tasks.size() > 1)
        for (auto& task : _factor_tasks)
        {
            task.taskFactor.reset(new LucasVMulFast(true));
            task.taskCheck.reset(new LucasVMulFast(true));
        }
    for (i = 0; i < input.factors().size(); i++)
    {
        Giant& b = input.factors()[i].first;
        if (b == 2 && _negQ)
            continue;
        n = input.factors()[i].second;
        auto factor = std::find_if(_factor_tasks.begin(), _factor_tasks.end(), [&](auto& a) { return a.index == i; });
        if (factor != _factor_tasks.end())
            n--;

        if (b.bitlen() < 32)
        {
            int prime = (int)b.data()[0];
            int index = 0;
            if (taskV && n > 0)
                index = taskV->mul_prime(prime, n);
            if (factor != _factor_tasks.end())
            {
                index = _taskCheck->mul_prime(prime, 1, index);
                if (_factor_tasks.size() > 1)
                    for (auto& task : _factor_tasks)
                    {
                        if (task.index == i)
                            task.taskCheck->mul_prime(prime, 1, index);
                        else
                            task.taskFactor->mul_prime(prime, 1, index);
                    }
            }
        }
        else
        {
            if (taskV && n > 0)
                taskV->mul_giant(b, n);
            if (factor != _factor_tasks.end())
            {
                _taskCheck->mul_giant(b, 1);
                if (_factor_tasks.size() > 1)
                    for (auto& task : _factor_tasks)
                    {
                        if (task.index == i)
                            task.taskCheck->mul_giant(b, 1);
                        else
                            task.taskFactor->mul_giant(b, 1);
                    }
            }
        }
    }

    for (_P = (_negQ ? 1 : 3); kronecker(_P*_P - (_negQ ? -4 : 4), exp) == 1; _P++);
    logging.progress().add_stage(_task->cost());
}

void Morrison::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging)
{
    if (!_task)
        return;

    if (logging.progress().param_int("P") != 0)
        _P = logging.progress().param_int("P");
    std::string sP = std::to_string(_P);

    File* checkpoint = nullptr;
    File* recoverypoint = nullptr;

    _success = false;
    _prime = false;
    bool restart = false;
    while (!_prime)
    {
        if (restart)
        {
            for (_P++; kronecker(_P*_P - (_negQ ? -4 : 4), *gwstate.N) == 1; _P++);
            sP = std::to_string(_P);

            logging.report_param("P", sP);
            logging.progress().add_stage(_task->cost());
            logging.progress().update(0, 0);
            logging.progress_save();

            if (checkpoint)
                checkpoint->clear();
            if (recoverypoint)
                recoverypoint->clear();
        }

        Giant tmp;
        tmp = _P*4;
        tmp *= _P*_P - (_negQ ? -4 : 4);
        tmp.gcd(*gwstate.N);
        if (tmp != 1)
        {
            _res64 = tmp.to_res64();
            logging.set_prefix("");
            logging.result(_prime, "%s is not prime. Factor RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), logging.progress().time_total());
            logging.result_save(input.input_text() + " is not prime. Factor RES64: " + _res64 + ", time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
            break;
        }

        logging.set_prefix("");
        logging.info("%sMorrison%s test of %s, P = %d, Q = %d, complexity = %d.\n", restart ? "Restarting " : "", _LLR ? " (LLR)" : "", input.display_text().data(), _P, _negQ ? -1 : 1, (int)logging.progress().cost_total());
        logging.set_prefix(input.display_text() + " ");
        if (gwstate.information_only)
            throw TaskAbortException();
        restart = true;

        checkpoint = file_checkpoint.add_child(sP, File::unique_fingerprint(file_checkpoint.fingerprint(), sP));
        if (dynamic_cast<LucasVMulFast*>(_task.get()) != nullptr)
            dynamic_cast<LucasVMulFast*>(_task.get())->init(&input, &gwstate, checkpoint, &logging, _P, _negQ);
        if (dynamic_cast<LucasUVMulFast*>(_task.get()) != nullptr)
        {
            recoverypoint = file_recoverypoint.add_child(sP, File::unique_fingerprint(file_recoverypoint.fingerprint(), sP));
            dynamic_cast<LucasUVMulFast*>(_task.get())->init(&input, &gwstate, checkpoint, recoverypoint, &logging, _P, _negQ);
        }
        _task->run();

        Giant* result = _task->result();
        bool result_parity = _task->result_parity();
        if (_taskCheck)
        {
            _taskCheck->init(&input, &gwstate, nullptr, &logging, _negQ);
            _taskCheck->init_state(new LucasVMulFast::State(0, 0, *result, result_parity));
            _taskCheck->run();
            result = _taskCheck->result();
        }
        logging.progress().next_stage();

        if (*result != (_negQ ? 0 : 2))
        {
            _res64 = result->to_res64();
            logging.set_prefix("");
            logging.result(_success, "%s is not a probable prime. Have you run Fermat test first? RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), logging.progress().time_total());
            logging.result_save(input.input_text() + " is not a probable prime. Have you run Fermat test first? RES64: " + _res64 + ", time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
            break;
        }
        _success = true;

        if (_factor_tasks.size() > 0)
        {
            Giant G, done;
            done = 1;
            std::string factors;
            if (_negQ)
                factors = "2";

            if (_factor_tasks.size() == 1)
            {
                G = std::move(*_task->result());
                if (!_negQ)
                    G -= 2;
                if (G == 0)
                    continue;
                factors += (!factors.empty() ? ", " : "") + input.factors()[_factor_tasks[0].index].first.to_string();
            }
            else
            {
                std::vector<Giant> Gs;
                for (auto& ftask : _factor_tasks)
                {
                    ftask.taskFactor->init(&input, &gwstate, nullptr, &logging, _negQ);
                    ftask.taskFactor->init_state(new LucasVMulFast::State(0, 0, *_task->result(), _task->result_parity()));
                    ftask.taskFactor->run();
                    if (ftask.taskFactor->state()->V() == (_negQ ? 0 : 2))
                    {
                        if (_all_factors)
                        {
                            done = 0;
                            break;
                        }
                        continue;
                    }
                    ftask.taskCheck->init(&input, &gwstate, nullptr, &logging, _negQ);
                    ftask.taskCheck->init_state(new LucasVMulFast::State(0, 0, *ftask.taskFactor->result(), ftask.taskFactor->result_parity()));
                    ftask.taskCheck->run();
                    if (ftask.taskCheck->state()->V() != (_negQ ? 0 : 2))
                    {
                        logging.warning("Arithmetic error, restarting.");
                        done = 0;
                        break;
                    }
                    Gs.push_back(std::move(*ftask.taskFactor->result()));
                    if (!_negQ)
                        Gs.back() -= 2;
                    done *= power(input.factors()[ftask.index].first, input.factors()[ftask.index].second);
                    factors += (!factors.empty() ? ", " : "") + input.factors()[ftask.index].first.to_string();
                }
                if (done == 0)
                    continue;
                int n = _negQ ? input.factors()[0].second : 0;
                if ((done.bitlen() + n)*2 + 10 < gwstate.N->bitlen())
                    continue;
                if ((done.bitlen() + n)*2 < gwstate.N->bitlen() + 10)
                {
                    done <<= n;
                    if (square(done) < *gwstate.N)
                        continue;
                }

                if (Gs.size() > 1)
                {
                    Product taskP(Gs.begin(), Gs.end());
                    taskP.init(&input, &gwstate, &logging);
                    taskP.run();
                    G = std::move(taskP.result());
                }
                else
                    G = std::move(Gs[0]);
            }

            logging.info("Checking gcd with factors {%s}.\n", factors.data());
            G.gcd(*gwstate.N);
            logging.progress().update(0, 0);
            if (G != 1) // Q=1: 19*2130-1, Q=-1: 225*5516-1
            {
                _res64 = G.to_res64();
                logging.set_prefix("");
                logging.result(_prime, "%s is not prime. Factor RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), logging.progress().time_total());
                logging.result_save(input.input_text() + " is not prime. Factor RES64: " + _res64 + ", time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
                break;
            }
        }
        _prime = true;
    }

    if (_prime)
    {
        logging.set_prefix("");
        logging.result(_prime, "%s is prime! Time: %.1f s.\n", input.display_text().data(), logging.progress().time_total());
        logging.result_save(input.input_text() + " is prime! Time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
    }

    if (checkpoint)
        checkpoint->clear();
    if (recoverypoint)
        recoverypoint->clear();
}

MorrisonGeneric::MorrisonGeneric(InputNum& input, Options& options, Logging& logging) : _options(options)
{
    bool CheckStrong = options.CheckStrong ? options.CheckStrong.value() : true;
    if (!CheckStrong)
        logging.error("The implementation of the Morrison test for such numbers requires the strong check.\n");
    CheckStrong = true;
    options.maxmulbyconst = 2;
    if (options.AllFactors)
        _all_factors = options.AllFactors.value();

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

    std::vector<std::unique_ptr<FactorTree>> factors;
    _dac_index.reserve(input.factors().size());
    PrimeIterator primes = PrimeIterator::get();
    for (int i = 0; i < input.factors().size(); i++)
    {
        auto& f = input.factors()[i];
        if (_done_factors.count(i) != 0)
        {
            tmp = f.first;
            tmp.power(f.second);
            tmp_done *= tmp;
        }
        if (f.first.bitlen() < 32)
        {
            int& prime = (int&)*f.first.data();
            while (primes.pos() < precomputed_DAC_S_d_len && *primes < prime)
                primes++;
            if (primes.pos() < precomputed_DAC_S_d_len)
            {
                GWASSERT(*primes == prime);
                GWASSERT(i == _dac_index.size());
                _dac_index.push_back((int)primes.pos());
            }
            else
            {
                int len = 60;
                _dac_index.push_back(-get_DAC_S_d(prime, (int)(prime/1.618) - 100, (int)(prime/1.618) + 100, &len));
            }
        }

        int power = f.second - 1;
        if (f.first == 2)
        {
            _negQ = (power > 0);
            if (!_negQ)
                factors.emplace_back(new FactorTree(f.first, i));
        }
        else
            factors.emplace_back(new FactorTree(f.first, i));
        if (power > 0)
        {
            tmp = f.first;
            tmp.power(power);
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

    _tree.reset(new FactorTree(factors));
    _tree->exp() = std::move(exp);

    _logging.reset(new SubLogging(logging, logging.level() + 1));
    _logging->progress().set_parent(nullptr);
    logging.progress().add_stage(input.bitlen()*std::log2(input.factors().size()));
}

void MorrisonGeneric::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging)
{
    if (logging.progress().param_int("P") != 0)
        _P = logging.progress().param_int("P");
    else
        _P = (_negQ ? 1 : 5);
    for (; kronecker(_P*_P - (_negQ ? -4 : 4), *gwstate.N) == 1; _P++);
    logging.info("Morrison test of %s, P = %d, Q = %d.\n", input.display_text().data(), _P, _negQ ? -1 : 1);
    if (gwstate.information_only)
        throw TaskAbortException();

    bool CheckStrong = _options.CheckStrong ? _options.CheckStrong.value() : true;
    CheckStrong = true;
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
        logging.progress().update((_all_factors ? 1 : 2)*_done.bitlen()/pct_bitlen/1000.0, 0);
        logging.heartbeat();
        last_progress = 0;
    };

    while (_res64.empty() && !_prime)
    {
        _success = false;
        _res64 = "";
        logging.set_prefix(input.display_text() + " ");
        if (_done > 1)
            logging.info("Restarting with %.1f%% of factors tested.\n", std::floor(_done.bitlen()/pct_bitlen)/10.0);

        Giant tmp;
        tmp = _P*4;
        tmp *= _P*_P - (_negQ ? -4 : 4);
        tmp.gcd(*gwstate.N);
        if (tmp != 1)
        {
            _res64 = tmp.to_res64();
            logging.set_prefix("");
            logging.result(_prime, "%s is not prime. Factor RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), logging.progress().time_total());
            logging.result_save(input.input_text() + " is not prime. Factor RES64: " + _res64 + ", time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
            break;
        }

        uint32_t fingerprint = File::unique_fingerprint(file_checkpoint.fingerprint(), std::to_string(input.factors().size()) + "." + std::to_string(_P));
        double last_write = 0;
        std::vector<FactorTree*> stack;
        std::vector<LucasVMulFast::State> stack_value;
        int task_num = 0;
        std::vector<int> factors;
        std::string factors_str;
        std::vector<Giant> G;
        auto test_G = [&]() {
            if (!G.empty())
            {
                if (G.size() > 1)
                {
                    bool a = Task::abort_flag();
                    Task::abort_reset();
                    Product taskP(G.begin(), G.end());
                    taskP.init(&input, &gwstate, _logging.get());
                    taskP.run();
                    tmp = std::move(taskP.result());
                    if (a)
                        Task::abort();
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
            }

            for (auto it = factors.begin(); it != factors.end(); it++)
            {
                logging.progress().param("factors") += "," + std::to_string(*it);
                _done_factors.insert(*it);
            }
            factors.clear();

            if (tmp_done != 1)
                _done *= tmp_done;
            tmp_done = 1;

            if ((!_all_factors || stack.empty()) && (_done.bitlen()*2 + 10 > gwstate.N->bitlen() &&
                (_done.bitlen()*2 > gwstate.N->bitlen() + 10 || square(_done) > *gwstate.N)))
            {
                _prime = true;
                return true;
            }
            return false;
        };
        std::string stack_file_prefix = std::to_string(_P) + ".s";
        auto write_values = [&]() {
            if (!G.empty())
            {
                if (test_G())
                    return true;
                logging.progress_save();
            }
            logging.debug("saving state to disk.\n");
            int i;
            for (i = 0; true; i++)
            {
                if (i < stack_value.size() && stack_value[i].is_written())
                    continue;
                std::string filename = file_checkpoint.filename() + "." + stack_file_prefix + std::to_string(i);
                auto it = file_checkpoint.children().begin();
                for (; it != file_checkpoint.children().end() && (*it)->filename() != filename; it++);
                if (it == file_checkpoint.children().end())
                {
                    if (i < stack_value.size())
                        continue;
                    else if (i == stack_value.size())
                    {
                        file_checkpoint.add_child(stack_file_prefix + std::to_string(i), fingerprint);
                        it = file_checkpoint.children().end() - 1;
                    }
                    else
                        break;
                }
                (*it)->clear();
                file_checkpoint.children().erase(it);
            }
            for (i = 0; i < stack_value.size(); i++)
            {
                if (stack_value[i].is_written())
                    continue;
                File* checkpoint = file_checkpoint.add_child(stack_file_prefix + std::to_string(i), fingerprint);
                checkpoint->write(stack_value[i]);
            }
            last_write = 0;
            return false;
        };

        stack.push_back(_tree.get());
        if (stack.back()->exp() == 1)
            task_num = 1;
        for (int i = 0; true; i++)
        {
            File* checkpoint = file_checkpoint.add_child(stack_file_prefix + std::to_string(i), fingerprint);
            stack_value.emplace_back();
            if (!checkpoint->read(stack_value.back()))
            {
                file_checkpoint.children().pop_back();
                stack_value.pop_back();
                break;
            }
            if (task_num > 0)
            {
                if (stack_value.back().index() == task_num + 1)
                    stack.push_back(stack.back()->left());
                else
                {
                    stack.back()->left()->exp().arithmetic().free(stack.back()->left()->exp());
                    stack.push_back(stack.back()->right());
                }
            }
            stack.back()->exp().arithmetic().free(stack.back()->exp());
            task_num = stack_value.back().index();
        }
        if (task_num > 0)
            stack.push_back(stack.back()->left());

        while (!stack.empty())
        {
            GWASSERT(stack_value.empty() || stack.back()->exp() != 1);
            if (!stack.back()->exp().empty() && stack.back()->exp() != 1)
            {
                if (stack.back()->index() >= 0)
                    GWASSERT(input.factors()[stack.back()->index()].first == stack.back()->exp());
                task_num++;
                int bitlen = stack.back()->exp().bitlen();
                int checks = _options.StrongCount ? _options.StrongCount.value() : 16;
                checks = bitlen > checks*1000 ? checks : bitlen/1000 + 1;

                std::unique_ptr<LucasVMul> cur_task;
                if (bitlen < 32)
                {
                    LucasVMulFast* task;
                    cur_task.reset(task = new LucasVMulFast(true));
                    if (stack.back()->index() >= 0 && stack.back()->index() < _dac_index.size())
                        task->mul_prime((int)*stack.back()->exp().data(), 1, _dac_index[stack.back()->index()]);
                    else
                        task->mul_giant(std::move(stack.back()->exp()), 1);
                }
                else if (stack_value.empty())
                {
                    if (CheckStrong)
                        cur_task.reset(new LucasUVMulFast(std::move(stack.back()->exp()), checks));
                }
                else
                {
                    if (CheckStrong)
                        cur_task.reset(new LucasUVMul(std::move(stack.back()->exp()), checks));
                }
                stack.back()->exp().arithmetic().free(stack.back()->exp());
                cur_task->set_error_check(!_options.CheckNear || _options.CheckNear.value(), _options.Check && _options.Check.value());
                _logging->progress() = Progress();
                _logging->progress().add_stage(cur_task->cost());
                _logging->set_prefix("stage " + std::to_string(task_num) + " ");

                File* checkpoint = nullptr;
                File* recoverypoint = nullptr;
                if (bitlen >= 32)
                {
                    std::string file_num = std::to_string(_P) + "." + std::to_string(task_num);
                    checkpoint = file_checkpoint.add_child(file_num, fingerprint);
                    recoverypoint = file_recoverypoint.add_child(file_num, fingerprint);
                }
                if (!stack_value.empty())
                {
                    if (LucasVMulFast* task = dynamic_cast<LucasVMulFast*>(cur_task.get()))
                        task->init(&input, &gwstate, checkpoint, _logging.get(), std::move(stack_value.back().V()), stack_value.back().parity());
                    if (LucasUVMul* task = dynamic_cast<LucasUVMul*>(cur_task.get()))
                        task->init(&input, &gwstate, checkpoint, recoverypoint, _logging.get(), std::move(stack_value.back().V()), stack_value.back().parity());
                }
                else
                {
                    if (LucasVMulFast* task = dynamic_cast<LucasVMulFast*>(cur_task.get()))
                        task->init(&input, &gwstate, checkpoint, _logging.get(), _P, _negQ);
                    if (LucasUVMul* task = dynamic_cast<LucasUVMul*>(cur_task.get()))
                        task->init(&input, &gwstate, checkpoint, recoverypoint, _logging.get(), _P, _negQ);
                }

                try
                {
                    cur_task->run();
                }
                catch (const TaskAbortException&)
                {
                    if (Task::abort_flag())
                    {
                        if (!stack_value.empty())
                            stack_value.back().V() = std::move(cur_task->P());
                        if (write_values())
                            break;
                    }
                    throw;
                }
                if (!stack_value.empty())
                    stack_value.back().V() = std::move(cur_task->P());
                stack_value.emplace_back(0, task_num, std::move(*cur_task->result()), cur_task->negativeQ() && cur_task->result_parity());
                last_progress += cur_task->timer();
                last_write += cur_task->timer();
                _logging->progress().next_stage();
                cur_task.reset();
                if (checkpoint)
                {
                    checkpoint->clear();
                    file_checkpoint.children().pop_back();
                }
                if (recoverypoint)
                {
                    recoverypoint->clear();
                    file_recoverypoint.children().pop_back();
                }

                if (stack.back()->is_factor())
                {
                    if (stack_value.back().V() != (_negQ ? 0 : 2) && !_success)
                    {
                        logging.set_prefix("");
                        _res64 = stack_value.back().V().to_res64();
                        logging.result(_success, "%s is not a probable prime. Have you run Fermat test first? RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), logging.progress().time_total());
                        logging.result_save(input.input_text() + " is not a probable prime. Have you run Fermat test first? RES64: " + _res64 + ", time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
                        break;
                    }
                    if (stack_value.back().V() != (_negQ ? 0 : 2))
                    {
                        logging.warning("Arithmetic error, restarting.");
                        break;
                    }
                    if (!_success && _negQ && _done_factors.count(0) == 0)
                    {
                        _done <<= input.factors()[0].second;
                        factors.push_back(0);
                    }
                    _success = true;
                    int index = stack.back()->index();
                    stack.pop_back();
                    stack_value.pop_back();
                    if (_done_factors.count(index) != 0)
                    {
                        _logging->info("Factor #%d has been tested already.\n", index);
                    }
                    else if (stack_value.back().V() != (_negQ ? 0 : 2))
                    {
                        _logging->info("Factor #%d added to gcd.\n", index);
                        G.push_back(std::move(stack_value.back().V()));
                        if (!_negQ)
                            G.back() -= 2;
                        factors.push_back(index);

                        auto& f = input.factors()[index];
                        tmp_done *= power(f.first, f.second);
                        if (tmp_done.size() > 8192)
                        {
                            _done *= tmp_done;
                            tmp_done = 1;
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
                        _logging->info("Factor #%d can't be tested with P=%d.\n", index, _P);
                    if (last_progress > Task::PROGRESS_TIME)
                        progress_factors();
                }
                else
                {
                    if (last_write > Task::DISK_WRITE_TIME)
                    {
                        if (write_values())
                            break;
                    }
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
        if (test_G())
            break;

        for (_P++; kronecker(_P*_P - (_negQ ? -4 : 4), *gwstate.N) == 1; _P++);
        logging.report_param("P", _P);
        logging.report_param("factors", "");
        _done = 1;
        _done_factors.clear();

        logging.progress_save();
        file_checkpoint.clear(true);
        file_checkpoint.children().clear();
        file_recoverypoint.clear(true);
        file_recoverypoint.children().clear();

        logging.set_prefix("");
        logging.info("Restarting Morrison test of %s, P = %d, Q = %d.\n", input.display_text().data(), _P, _negQ ? -1 : 1);
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
