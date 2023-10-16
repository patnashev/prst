
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
// For Q=1  gcd(U_{(N+1)/q}, N) = gcd(V_{2(N+1)/q} - 2, N)
// For Q=-1 gcd(U_{(N+1)/q}, N) = gcd(V_{(N+1)/2q}, N) due to BLS proof of Theorem 14.
// For Q=-1 factor 2 is tested for free.

Morrison::Morrison(InputNum& input, Params& params, Logging& logging)
{
    int i;
    Giant tmp;
    std::vector<std::pair<Giant, int>> factors;

    bool CheckStrong = params.CheckStrong ? params.CheckStrong.value() : false;
    if (params.AllFactors)
        _all_factors = params.AllFactors.value();

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
            n = factor.second;
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
    if (n > 1)
    {
        _negQ = true;
        n--;
    }
    else
    {
        _negQ = false;
        n++;
    }

    if (CheckStrong)
    {
        exp_morrison <<= n;
        int checks = params.StrongCount ? params.StrongCount.value() : 16;
        _task.reset(new LucasUVMulFast(std::move(exp_morrison), checks));
        params.maxmulbyconst = 2;
    }
    else
        taskV->mul_prime(2, n);

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
        if (b == 2)
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
                        if (task.index == i)
                            task.taskCheck->mul_prime(prime, 1, index);
                        else
                            task.taskFactor->mul_prime(prime, 1, index);
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
                        if (task.index == i)
                            task.taskCheck->mul_giant(b, 1);
                        else
                            task.taskFactor->mul_giant(b, 1);
            }
        }
    }

    for (_P = (_negQ ? 1 : 3); kronecker(_P*_P - (_negQ ? -4 : 4), exp) == 1; _P++);
    logging.progress().add_stage(_task->cost());
}

void Morrison::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, File& file_params, Logging& logging)
{
    if (!_task)
        return;

    File* checkpoint = nullptr;
    File* recoverypoint = nullptr;
    std::unique_ptr<Reader> reader(file_params.get_reader());
    if (reader && reader->type() == 7)
        reader->read(_P);

    _success = false;
    _prime = false;
    bool restart = false;
    while (!_prime)
    {
        if (restart)
        {
            for (_P++; kronecker(_P*_P - (_negQ ? -4 : 4), *gwstate.N) == 1; _P++);
            logging.progress().add_stage(_task->cost());

            std::unique_ptr<Writer> writer(file_params.get_writer(7, 0));
            writer->write(_P);
            file_params.commit_writer(*writer);
            if (checkpoint)
                checkpoint->clear();
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
            exit(0);
        restart = true;

        checkpoint = file_checkpoint.add_child(std::to_string(_P), File::unique_fingerprint(file_checkpoint.fingerprint(), std::to_string(_P)));
        if (dynamic_cast<LucasVMulFast*>(_task.get()) != nullptr)
            dynamic_cast<LucasVMulFast*>(_task.get())->init(&input, &gwstate, checkpoint, &logging, _P, _negQ);
        if (dynamic_cast<LucasUVMulFast*>(_task.get()) != nullptr)
        {
            recoverypoint = file_recoverypoint.add_child(std::to_string(_P), File::unique_fingerprint(file_recoverypoint.fingerprint(), std::to_string(_P)));
            dynamic_cast<LucasUVMulFast*>(_task.get())->init(&input, &gwstate, checkpoint, recoverypoint, &logging, _P, _negQ);
        }
        _task->run();

        Giant* result = _task->result();
        if (_taskCheck)
        {
            _taskCheck->init(&input, &gwstate, nullptr, &logging, _negQ);
            _taskCheck->init_state(new LucasVMulFast::State(0, 0, *result, false));
            _taskCheck->run();
            result = &_taskCheck->state()->V();
        }
        logging.progress().next_stage();

        if (*result != (_negQ ? 0 : 2))
        {
            LucasVMulFast taskPRP(true);
            if (_negQ)
            {
                _res64 = result->to_res64();
                taskPRP.mul_prime(2, 2, 0);
                taskPRP.init(&input, &gwstate, nullptr, &logging, _negQ);
                taskPRP.init_state(new LucasVMulFast::State(0, 0, *result, false));
                taskPRP.run();
                result = &taskPRP.state()->V();
            }
            if (*result != 2)
            {
                _res64 = result->to_res64();
                logging.set_prefix("");
                logging.result(_success, "%s is not a probable prime. Have you run Fermat test first? RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), logging.progress().time_total());
                logging.result_save(input.input_text() + " is not a probable prime. Have you run Fermat test first? RES64: " + _res64 + ", time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
                break;
            }
            else
            {
                _success = true;
                logging.set_prefix("");
                logging.result(_prime, "%s is not prime. RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), logging.progress().time_total());
                logging.result_save(input.input_text() + " is not prime. RES64: " + _res64 + ", time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
                break;
            }
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
                    ftask.taskFactor->init_state(new LucasVMulFast::State(0, 0, *_task->result(), false));
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
                    ftask.taskCheck->init_state(new LucasVMulFast::State(0, 0, ftask.taskFactor->state()->V(), ftask.taskFactor->state()->parity()));
                    ftask.taskCheck->run();
                    if (ftask.taskCheck->state()->V() != (_negQ ? 0 : 2))
                    {
                        logging.warning("Arithmetic error, restarting.");
                        done = 0;
                        break;
                    }
                    Gs.push_back(std::move(ftask.taskFactor->state()->V()));
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

            logging.result(false, "Checking gcd with factors {%s}.\n", factors.data());
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
    file_params.clear();
}
