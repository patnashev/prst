
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

    if (params.CheckStrong)
        logging.warning("Strong check is not implemented in Morrison test. Use -fermat first.\n");

    tmp = 1;
    factors.reserve(input.b_factors().size());
    for (i = 0; i < input.b_factors().size(); i++)
        if (input.b_factors()[i].first != 2)
            factors.emplace_back(power(input.b_factors()[i].first, input.b_factors()[i].second), i);
        else
        {
            tmp = power(input.b_factors()[i].first, input.b_factors()[i].second);
            _factors = "2";
        }
    std::sort(factors.begin(), factors.end(), [](const std::pair<Giant, int>& a, const std::pair<Giant, int>& b) { return a.first > b.first; });
    for (i = 0; i < factors.size() && tmp*tmp < input.gb(); i++)
    {
        tmp *= factors[i].first;
        _factors += (!_factors.empty() ? ", " : "") + input.b_factors()[factors[i].second].first.to_string();
    }

    Giant N = input.value();
    _negQ = N.bit(0) && N.bit(1);
    _task.reset(new LucasMul(_negQ));
    if (!_negQ)
        _task->mul_prime(2, 1, 0);
    if (_factor_tasks.size() > 0)
        _taskCheck.reset(new LucasMul(_negQ));

    if (_factor_tasks.size() > 1)
        for (i = 0; i < _factor_tasks.size(); i++)
        {
            _factor_tasks[i].taskFactor.reset(new LucasMul(_negQ));
            _factor_tasks[i].taskCheck.reset(new LucasMul(_negQ));
        }
    
    bool div2 = false;
    if (input.gk() != 1)
        if (!input.gk().bit(0) && _negQ)
        {
            div2 = true;
            tmp = input.gk();
            tmp >>= 1;
            _task->mul_giant(std::move(tmp), 1);
        }
        else
            _task->mul_giant(input.gk(), 1);

    
    for (i = 0; i < input.b_factors().size(); i++)
    {
        Giant& b = input.b_factors()[i].first;
        int n = input.n()*input.b_factors()[i].second;
        if (!div2 && _negQ && b == 2)
        {
            div2 = true;
            n--;
        }
        auto factor = std::find_if(_factor_tasks.begin(), _factor_tasks.end(), [&](auto& a) { return a.index == i; });
        if (factor != _factor_tasks.end())
            n--;

        if (b.bitlen() < 32)
        {
            int prime = (int)b.data()[0];
            int index = n > 0 ? _task->mul_prime(prime, n) : 0;
            if (factor != _factor_tasks.end())
            {
                _taskCheck->mul_prime(prime, 1, index);
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
            _task->mul_giant(b, n);
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
    
    for (_P = 3; kronecker(_P*_P - (_negQ ? -4 : 4), N) == 1; _P++);
    logging.progress().add_stage(_task->cost());

    GWASSERT(div2 == _negQ);
}

void Morrison::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_params, Logging& logging)
{
    File* checkpoint = nullptr;
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
        restart = true;

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
        logging.info("Morrison test of %s, P = %d, Q = %d, factors = {%s}, complexity = %d.\n", input.display_text().data(), _P, _negQ ? -1 : 1, _factors.data(), (int)logging.progress().cost_total());
        logging.set_prefix(input.display_text() + " ");
        if (gwstate.information_only)
            exit(0);

        checkpoint = file_checkpoint.add_child(std::to_string(_P), File::unique_fingerprint(file_checkpoint.fingerprint(), std::to_string(_P)));
        _task->init(&input, &gwstate, checkpoint, &logging, _P);
        _task->run();

        LucasMul::State* state = _task->state();
        if (_taskCheck)
        {
            _taskCheck->init(&input, &gwstate, nullptr, &logging);
            _taskCheck->init_state(new LucasMul::State(0, 0, _task->state()->V(), _task->state()->parity()));
            _taskCheck->run();
            state = _taskCheck->state();
        }
        logging.progress().next_stage();

        if (state->V() != (_negQ ? 0 : 2))
        {
            LucasMul taskPRP(_negQ);
            if (_negQ)
            {
                _res64 = state->V().to_res64();
                taskPRP.mul_prime(2, 2, 0);
                taskPRP.init(&input, &gwstate, nullptr, &logging);
                taskPRP.init_state(new LucasMul::State(0, 0, state->V(), state->parity()));
                taskPRP.run();
                state = taskPRP.state();
            }
            if (state->V() != 2)
            {
                _res64 = state->V().to_res64();
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
            Giant G;
            if (_factor_tasks.size() == 1)
            {
                G = std::move(_task->state()->V());
                if (!_negQ)
                    G -= 2;
                if (G == 0)
                    continue;
            }
            else
            {
                std::vector<Giant> Gs;
                for (auto& ftask : _factor_tasks)
                {
                    ftask.taskFactor->init(&input, &gwstate, nullptr, &logging);
                    ftask.taskFactor->init_state(new LucasMul::State(0, 0, _task->state()->V(), _task->state()->parity()));
                    ftask.taskFactor->run();
                    if (ftask.taskFactor->state()->V() == (_negQ ? 0 : 2))
                    {
                        Gs.clear();
                        break;
                    }
                    ftask.taskCheck->init(&input, &gwstate, nullptr, &logging);
                    ftask.taskCheck->init_state(new LucasMul::State(0, 0, ftask.taskFactor->state()->V(), ftask.taskFactor->state()->parity()));
                    ftask.taskCheck->run();
                    if (ftask.taskCheck->state()->V() != (_negQ ? 0 : 2))
                    {
                        logging.warning("Arithmetic error, restarting.");
                        Gs.clear();
                        break;
                    }
                    Gs.push_back(std::move(ftask.taskFactor->state()->V()));
                    if (!_negQ)
                        Gs.back() -= 2;
                }
                if (Gs.empty())
                    continue;

                Product taskP(Gs.begin(), Gs.end());
                taskP.init(&input, &gwstate, &logging);
                taskP.run();
                G = std::move(taskP.result());
            }

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

int LucasMul::mul_prime(int prime, int n, int index)
{
    if (prime < 14)
        index = 0;
    else if (index == 0)
    {
        if (prime < precomputed_DAC_S_d_len*(log(precomputed_DAC_S_d_len) + log(log(precomputed_DAC_S_d_len))))
            for (auto it = PrimeIterator::get(); *it != prime && index < precomputed_DAC_S_d_len; index++, it++);
        if (index == 0 || index >= precomputed_DAC_S_d_len)
        {
            int len = 60;
            index = -get_DAC_S_d(prime, (int)(prime/1.618) - 100, (int)(prime/1.618) + 100, &len);
        }
    }
    _primes.emplace_back(prime, n, index);
    _progress.reset();
    return index;
}

double LucasMul::cost()
{
    if (!_progress)
        progress_init();
    return _progress->cost_total();
}

double LucasMul::progress()
{
    if (!_progress)
        progress_init();
    if (_state)
        _progress->update(_state->iteration()/(double)iterations(), 0);
    return _progress->progress_total();
}

void LucasMul::progress_init()
{
    _progress.reset(new Progress());
    for (auto& giant : _giants)
        _progress->add_stage(2*giant.first.bitlen()*giant.second);
    for (auto& prime : _primes)
    {
        int len = 60;
        switch (std::get<0>(prime))
        {
        case 2:
            len = 1; break;
        case 3:
            len = 2; break;
        case 5:
            len = 3; break;
        case 7:
            len = 4; break;
        case 11:
            len = 5; break;
        case 13:
            len = 6; break;
        default:
            int d = std::get<2>(prime) < 0 ? -std::get<2>(prime) : precomputed_DAC_S_d[std::get<2>(prime)];
            get_DAC_S_d(std::get<0>(prime), d, d + 1, &len);
        }
        _progress->add_stage(len*std::get<1>(prime));
    }
}

void LucasMul::init(InputNum* input, GWState* gwstate, File* file, Logging* logging)
{
    InputTask::init(input, gwstate, file, nullptr, logging, 0);
    if (_error_check)
        _logging->info("max roundoff check enabled.\n");
    State* state = read_state<State>(file);
    if (state != nullptr)
        init_state(state);
}

void LucasMul::init_state(State* state)
{
    if (state == nullptr)
    {
        _state.reset();
        return;
    }
    _state.reset(state);
    if (!_progress)
        progress_init();
    if (_progress->cur_stage() > state->index())
        _progress->reset();
    while (_progress->cur_stage() < state->index())
        _progress->next_stage();
    if (state->index() < _giants.size())
        _iterations = std::get<1>(_giants[state->index()]);
    else
        _iterations = std::get<1>(_primes[state->index() - _giants.size()]);
    _logging->progress().update(progress(), (int)_gwstate->handle.fft_count/2);
    if (state->index() > 0 || _state->iteration() > 0)
        _logging->info("restarting at %.1f%%.\n", 100.0*progress());
}

void LucasMul::setup()
{

}

void LucasMul::release()
{

}

void LucasMul::execute()
{
    int i;
    LucasVArithmetic lucas(gw(), _negativeQ);

    GWASSERT(state() != nullptr);
    LucasV Vn(lucas, state()->V(), state()->parity());
    LucasV Vn1(lucas);

    int index = state()->index();
    i = state()->iteration();

    while (index < _giants.size())
    {
        _iterations = _giants[index].second;
        _state_update_period = MULS_PER_STATE_UPDATE/(2*_giants[index].first.bitlen());
        for (; i < _iterations; i++, commit_execute<State>(i, index, Vn.V(), Vn.parity()))
            lucas.mul(Vn, _giants[index].first, Vn, Vn1);

        i = 0;
        index++;
        _progress->next_stage();
        set_state<State>(0, index, std::move(state()->V()), state()->parity());
    }

    while (index - _giants.size() < _primes.size())
    {
        auto& prime = _primes[index - _giants.size()];
        _iterations = std::get<1>(prime);
        _state_update_period = MULS_PER_STATE_UPDATE/(int)(_progress->costs()[index]/_iterations);
        for (; i < _iterations; i++, commit_execute<State>(i, index, Vn.V(), Vn.parity()))
            lucas.mul(Vn, std::get<0>(prime), std::get<2>(prime), Vn);

        i = 0;
        index++;
        _progress->next_stage();
        set_state<State>(0, index, std::move(state()->V()), state()->parity());
    }
}
