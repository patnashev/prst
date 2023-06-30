
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
    if (type() != POCKLINGTON)
        return;

    _tasks.reserve(input.factors().size());
    for (auto i = 0; i < input.factors().size(); i++)
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

    if (params.AllFactors)
        _all_factors = params.AllFactors.value();
}

void Pocklington::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof)
{
    if (type() != POCKLINGTON)
    {
        Fermat::run(input, gwstate, file_checkpoint, file_recoverypoint, logging, proof);
        return;
    }

    logging.info("Pocklington test of %s, a = %d, complexity = %d.\n", input.display_text().data(), _a, (int)logging.progress().cost_total());
    Fermat::run(input, gwstate, file_checkpoint, file_recoverypoint, logging, proof);

    Giant done;
    done = 1;

    File* checkpoint = &file_checkpoint;
    File* recoverypoint = &file_recoverypoint;
    Giant tmp;
    while (!_tasks.empty())
    {
        if (!success())
            return;
        std::vector<Giant> G;
        G.reserve(_tasks.size());
        std::string factors;

        for (auto it = _tasks.begin(); it != _tasks.end(); )
        {
            if (it->taskFactor)
            {
                it->taskFactor->init(&input, &gwstate, &logging, std::move(*_task->result()));
                it->taskFactor->run();
                *_task->result() = std::move(it->taskFactor->X0());
                tmp = std::move(*it->taskFactor->result());

                it->taskCheck->init(&input, &gwstate, &logging, std::move(tmp));
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
                done *= power(input.factors()[it->index].first, input.factors()[it->index].second);
                factors += (!factors.empty() ? ", " : "") + input.factors()[it->index].first.to_string();
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
            logging.info("Checking gcd with factors {%s}.\n", factors.data());
            logging.set_prefix("");
            tmp.gcd(*gwstate.N);
            if (tmp != 1)
            {
                _res64 = tmp.to_res64();
                logging.result(_prime, "%s is not prime. Factor RES64: %s.\n", input.display_text().data(), _res64.data());
                logging.result_save(input.input_text() + " is not prime. Factor RES64: " + _res64 + ".\n");
                break;
            }
        }

        if (_tasks.empty() || (!_all_factors && done*done >= *gwstate.N))
        {
            _prime = true;
            break;
        }
        else
        {
            if (proof == nullptr)
            {
                logging.progress().add_stage(_task->cost());
                logging.progress().update(0, (int)gwstate.handle.fft_count/2);
                logging.progress_save();
                if (dynamic_cast<StrongCheckMultipointExp*>(_task.get()) != nullptr)
                    recoverypoint->write(*_task->state());
                else
                    checkpoint->write(*_task->state());
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

            logging.warning("Restarting Pocklington test of %s, a = %d.\n", input.display_text().data(), _a);
            checkpoint = file_checkpoint.add_child(sa, File::unique_fingerprint(file_checkpoint.fingerprint(), sa));
            recoverypoint = file_recoverypoint.add_child(sa, File::unique_fingerprint(file_recoverypoint.fingerprint(), sa));
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
