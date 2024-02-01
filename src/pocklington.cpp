
#include <algorithm>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "pocklington.h"
#include "integer.h"

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
        std::string factors;

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
                logging.report_param("factor" + std::to_string(it->index), "done");
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
            logging.result(false, "Checking gcd with factors {%s}.\n", factors.data());
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

            logging.warning("Restarting Pocklington test of %s, a = %d.\n", input.display_text().data(), _a);
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
