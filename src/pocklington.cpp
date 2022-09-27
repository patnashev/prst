
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
    for (i = 0; i < j; i++)
        _tasks.emplace_back(SlowExp(input.gb()/input.b_factors()[factors[i].second].first), factors[i].second);

    if (input.is_base2())
    {
        InputNum input_k;
        InputNum input_base2;
        input.to_base2(input_k, input_base2);
        _a = genProthBase(input_base2.gk(), input_base2.n());
        if (dynamic_cast<FastExp*>(_task.get()) != nullptr)
            params.maxmulbyconst = _a;
    }
}

void Pocklington::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof)
{
    std::string factors;
    for (auto it = _tasks.begin(); it != _tasks.end(); it++)
        factors += (!factors.empty() ? ", " : "") + input.b_factors()[it->second].first.to_string();
    logging.info("Pocklington test of %s, a = %d, factors = {%s}.\n", input.display_text().data(), _a, factors.data());
    Fermat::run(input, gwstate, file_checkpoint, file_recoverypoint, logging, proof);

    File* checkpoint = &file_checkpoint;
    File* recoverypoint = &file_recoverypoint;
    Giant tmp;
    int total = (int)_tasks.size();
    double timer = 0;
    std::vector<Giant> G;
    G.reserve(_tasks.size());
    while (_tasks.size() > 0)
    {
        if (!success())
            return;
        for (auto it = _tasks.begin(); it != _tasks.end(); )
        {
            SlowExp& task = it->first;
            if (task.exp() == 1)
                tmp = _Xm1;
            else
            {
                task.init(&input, &gwstate, nullptr, &logging, std::move(_Xm1));
                task.run();
                timer += task.timer();
                _Xm1 = std::move(task.X0());
                tmp = std::move(task.state()->X());
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

        if (_tasks.size() > 0)
        {
            if (proof == nullptr)
            {
                logging.progress().add_stage(_task->cost());
                logging.progress().update(0, (int)gwstate.handle.fft_count/2);
                logging.progress_save();
                BaseExp::State state(_task->state()->iteration(), std::move(_Xm1));
                if (dynamic_cast<GerbiczCheckMultipointExp*>(_task.get()) != nullptr)
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
            if (dynamic_cast<FastExp*>(_task.get()) != nullptr)
            {
                double fft_count = gwstate.handle.fft_count;
                gwstate.done();
                gwstate.maxmulbyconst = _a;
                input.setup(gwstate);
                gwstate.handle.fft_count = fft_count;
            }

            factors = "";
            for (auto it = _tasks.begin(); it != _tasks.end(); it++)
                factors += (!factors.empty() ? ", " : "") + input.b_factors()[it->second].first.to_string();
            logging.warning("Restarting Pocklington test of %s, a = %d, factors = {%s}.\n", input.display_text().data(), _a, factors.data());
            checkpoint = file_checkpoint.add_child(sa, File::unique_fingerprint(file_checkpoint.fingerprint(), sa));
            recoverypoint = file_recoverypoint.add_child(sa, File::unique_fingerprint(file_recoverypoint.fingerprint(), sa));
            Fermat::run(input, gwstate, *checkpoint, *recoverypoint, logging, nullptr);
        }
    }
    GWASSERT(G.size() == total);

    if (G.size() > 1)
    {
        Product taskP(G.begin(), G.end());
        taskP.init(&input, &gwstate, nullptr, &logging);
        taskP.run();
        timer += taskP.timer();
        tmp = std::move(taskP.state()->X());
    }
    else
        tmp = std::move(G[0]);
    tmp.gcd(*gwstate.N);
    if (tmp == 1)
    {
        _success = true;
        logging.result(_success, "%s is prime! Time: %.1f s.\n", input.display_text().data(), timer);
        logging.result_save(input.input_text() + " is prime! Time: " + std::to_string((int)timer) + " s.\n");
    }
    else
    {
        _success = false;
        _res64 = tmp.to_res64();
        logging.result(_success, "%s is not prime. Factor RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), timer);
        logging.result_save(input.input_text() + " is not prime. Factor RES64: " + _res64 + ", time: " + std::to_string((int)timer) + " s.\n");
    }

    file_checkpoint.clear(true);
    file_recoverypoint.clear(true);
}

template<class IT>
void Product<IT>::init(InputNum* input, GWState* gwstate, File* file, Logging* logging)
{
    BaseExp::init(input, gwstate, file, read_state<State>(file), logging, (int)(_last - _first));
    _state_update_period = MULS_PER_STATE_UPDATE;
}

template<class IT>
void Product<IT>::execute()
{
    GWNum P(gw());
    GWNum X(gw());
    P = *_first;
    int i = 1;
    for (IT it = _first + 1; it != _last; it++, i++, commit_execute<State>(i, P))
    {
        X = *it;
        gw().carefully().mul(X, P, P, 0);
    }

    done();
}
