
#include <algorithm>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "pocklington.h"
#include "integer.h"

using namespace arithmetic;

int genProthBase(Giant& k, uint32_t n);

Pocklington::Pocklington(InputNum& input, Params& params, Logging& logging, Proof* proof) : Fermat(Fermat::POCKLINGTON, input, params, logging, proof)
{
    int i, j;
    std::vector<std::pair<Giant, int>> factors;
    factors.reserve(input.b_factors().size());
    for (i = 0; i < input.b_factors().size(); i++)
        factors.emplace_back(power(input.b_factors()[i].first, input.b_factors()[i].second), i);
    std::sort(factors.begin(), factors.end(), [](const std::pair<Giant, int>& a, const std::pair<Giant, int>& b) { return a.first > b.first; });
    Giant tmp;
    tmp = 1;
    for (j = 0; j < factors.size() && tmp*tmp < input.gb(); j++)
        tmp *= factors[j].first;
    _tasks.reserve(j);
    for (i = 0; i < j; i++)
    {
        _tasks.emplace_back(factors[i].second);
        Giant& b = input.b_factors()[factors[i].second].first;
        _factors += (!_factors.empty() ? ", " : "") + b.to_string();
        tmp = input.gb()/b;
        if (tmp != 1)
        {
            _tasks.back().taskFactor.reset(new CarefulExp(std::move(tmp)));
            _tasks.back().taskCheck.reset(new CarefulExp(b));
        }
    }

    if (input.is_base2())
    {
        _input_k.reset(new InputNum());
        _input_base2.reset(new InputNum());
        input.to_base2(*_input_k, *_input_base2);
        _a = genProthBase(_input_base2->gk(), _input_base2->n());
        if (!_task->smooth())
            params.maxmulbyconst = _a;
    }
}

void Pocklington::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof)
{
    if (_input_base2)
        logging.info("Proth test of %s = %s, a = %d, complexity = %d.\n", input.display_text().data(), _input_base2->display_text().data(), _a, (int)logging.progress().cost_total());
    else
        logging.info("Pocklington test of %s, a = %d, factors = {%s}, complexity = %d.\n", input.display_text().data(), _a, _factors.data(), (int)logging.progress().cost_total());
    Fermat::run(input, gwstate, file_checkpoint, file_recoverypoint, logging, proof);

    File* checkpoint = &file_checkpoint;
    File* recoverypoint = &file_recoverypoint;
    Giant tmp;
    while (_tasks.size() > 0)
    {
        if (!success())
            return;
        std::vector<Giant> G;
        G.reserve(_tasks.size());

        for (auto it = _tasks.begin(); it != _tasks.end(); )
        {
            if (it->taskFactor)
            {
                it->taskFactor->init(&input, &gwstate, &logging, std::move(_Xm1));
                it->taskFactor->run();
                _Xm1 = std::move(it->taskFactor->X0());
                tmp = std::move(it->taskFactor->state()->X());

                it->taskCheck->init(&input, &gwstate, &logging, std::move(tmp));
                it->taskCheck->run();
                if (it->taskCheck->state()->X() != 1)
                {
                    logging.warning("Arithmetic error, restarting.");
                    continue;
                }
                tmp = std::move(it->taskCheck->X0());
            }
            else
                tmp = _Xm1;
            if (_input_base2)
            {
                tmp += 1;
                if (tmp != 0 && tmp != *gwstate.N)
                {
                    _res64 = tmp.to_res64();
                    logging.result(_prime, "%s is not prime. Proth RES64: %s.\n", input.display_text().data(), _res64.data());
                    logging.result_save(input.input_text() + " is not prime. Proth RES64: " + _res64 + ".\n");
                }
                break;
            }
            if (tmp != 1)
            {
                tmp -= 1;
                G.emplace_back(std::move(tmp));
                it = _tasks.erase(it);
            }
            else
                it++;
        }
        if (_input_base2)
            break;

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
            tmp.gcd(*gwstate.N);
            if (tmp != 1)
            {
                _res64 = tmp.to_res64();
                logging.result(_prime, "%s is not prime. Factor RES64: %s.\n", input.display_text().data(), _res64.data());
                logging.result_save(input.input_text() + " is not prime. Factor RES64: " + _res64 + ".\n");
                break;
            }
        }

        if (_tasks.size() > 0)
        {
            if (proof == nullptr)
            {
                logging.progress().add_stage(_task->cost());
                logging.progress().update(0, (int)gwstate.handle.fft_count/2);
                logging.progress_save();
                BaseExp::State state(_task->smooth() ? _n : _task->exp().bitlen() - 1, std::move(_Xm1));
                if (dynamic_cast<StrongCheckMultipointExp*>(_task.get()) != nullptr)
                    recoverypoint->write(state);
                else
                    checkpoint->write(state);
            }
            else
            {
                logging.error("Pocklington test needs to restart, disable proofs to proceed.\n", input.display_text().data(), _a);
                throw TaskAbortException();
            }

            PrimeIterator primes = PrimeIterator::get();
            for (; *primes <= _a; primes++);
            _a = *primes;
            std::string sa = std::to_string(_a);
            if (!_task->smooth())
            {
                double fft_count = gwstate.handle.fft_count;
                gwstate.done();
                gwstate.maxmulbyconst = _a;
                input.setup(gwstate);
                gwstate.handle.fft_count = fft_count;
            }

            std::string factors;
            for (auto it = _tasks.begin(); it != _tasks.end(); it++)
                factors += (!factors.empty() ? ", " : "") + input.b_factors()[it->index].first.to_string();
            logging.warning("Restarting Pocklington test of %s, a = %d, factors = {%s}.\n", input.display_text().data(), _a, factors.data());
            checkpoint = file_checkpoint.add_child(sa, File::unique_fingerprint(file_checkpoint.fingerprint(), sa));
            recoverypoint = file_recoverypoint.add_child(sa, File::unique_fingerprint(file_recoverypoint.fingerprint(), sa));
            Fermat::run(input, gwstate, *checkpoint, *recoverypoint, logging, nullptr);
        }
    }

    if (_success)
    {
        _prime = true;
        logging.result(_prime, "%s is prime!\n", input.display_text().data());
        logging.result_save(input.input_text() + " is prime!\n");
    }

    file_checkpoint.clear(true);
    file_recoverypoint.clear(true);
}
