
#include <cmath>
#include <algorithm>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "lucasmul.h"
#include "integer.h"

using namespace arithmetic;

void LucasVMul::done()
{
    InputTask::done();
    if (_gwstate->need_mod())
        _gwstate->mod(*result(), *result());
}

int LucasVMulFast::mul_prime(int prime, int n, int index)
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

double LucasVMulFast::cost()
{
    if (!_progress)
        progress_init();
    return _progress->cost_total();
}

double LucasVMulFast::progress()
{
    if (!_progress)
        progress_init();
    if (_state)
        _progress->update(_state->iteration()/(double)iterations(), 0);
    return _progress->progress_total();
}

void LucasVMulFast::progress_init()
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

void LucasVMulFast::init(InputNum* input, GWState* gwstate, File* file, Logging* logging)
{
    InputTask::init(input, gwstate, file, nullptr, logging, 0);
    if (_error_check && !_carefully)
        _logging->info("max roundoff check enabled.\n");
    State* state = read_state<State>(file);
    if (state != nullptr)
        init_state(state);
}

void LucasVMulFast::init_state(State* state)
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
    _logging->progress().update(progress(), 0);
    if (state->index() > 0 || _state->iteration() > 0)
        _logging->info("restarting at %.1f%%.\n", 100.0*progress());
    if (result() != nullptr && _gwstate->need_mod())
        _gwstate->mod(*result(), *result());
}

void LucasVMulFast::setup()
{

}

void LucasVMulFast::release()
{

}

void LucasVMulFast::execute()
{
    if (result() != nullptr)
        return;

    int i;
    LucasVArithmetic lucas(_carefully ? gw().carefully() : gw(), _negativeQ);

    GWASSERT(state() != nullptr);
    LucasV Vn(lucas, state()->V(), state()->parity());
    LucasV Vn1(lucas);

    int index = state()->index();
    i = state()->iteration();

    if (!_carefully && index == 0 && i == 0)
        gwset_carefully_count(gw().gwdata(), 30);

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

    done();
}

void LucasUVMulFast::Gerbicz_params(int iters, int& L, int &L2)
{
    int i;
    L = (int)std::sqrt(iters);
    L2 = iters + (L - iters%L)%L;
    L = L2/L;
    for (i = L + 1; i*i < 2*iters; i++)
        if (L2 > iters + (i - iters%i)%i)
        {
            L = i;
            L2 = iters + (i - iters%i)%i;
        }
}

void LucasUVMulFast::init(InputNum* input, GWState* gwstate, File* file, File* file_recovery, Logging* logging)
{
    InputTask::init(input, gwstate, file, read_state<StrongCheckState>(file), logging, _exp.bitlen());
    _state_update_period = MULS_PER_STATE_UPDATE/2;
    _logging->info("Gerbicz-Li check enabled, L2 = %d*%d.\n", _L, _L2/_L);
    _logging->report_param("L", _L);
    _logging->report_param("L2", _L2);
    if (_error_check)
        _logging->info("max roundoff check enabled.\n");
    _file_recovery = file_recovery;
    _state_recovery.reset();
    State* state_recovery = read_state<State>(file_recovery);
    if (state_recovery != nullptr)
        init_state(state_recovery);
}

void LucasUVMulFast::init_state(State* state)
{
    if (state == nullptr)
    {
        _state_recovery.reset();
        _state.reset();
        return;
    }
    _state_recovery.reset(state);
    _recovery_op = 0;
    if (!_state || (state_check() != nullptr ? state_check()->recovery() : _state->iteration()) != _state_recovery->iteration())
    {
        _state.reset(new TaskState(5));
        _state->set(_state_recovery->iteration());
    }
    _logging->progress().update(progress(), 0);
    if (_state->iteration() > 0)
        _logging->info("restarting at %.1f%%.\n", 100.0*_state->iteration()/iterations());
    if (result() != nullptr && _gwstate->need_mod())
        _gwstate->mod(*result(), *result());
}

void LucasUVMulFast::write_state()
{
    if (_file_recovery != nullptr && _state_recovery && !_state_recovery->is_written())
        _file_recovery->write(*_state_recovery);
    Task::write_state();
}

void LucasUVMulFast::setup()
{
    if (!_arithmetic)
        _arithmetic.reset(new LucasUVArithmetic(gw(), _P, _negativeQ));
    if (!_X)
        _X.reset(new LucasUV(arithmetic()));
    if (!_R)
        _R.reset(new LucasUV(arithmetic()));
    if (!_D)
        _D.reset(new LucasUV(arithmetic()));
    GWASSERT(arithmetic().max_small() > 0);
    for (_W = 2; (1 << _W) - 1 <= arithmetic().max_small(); _W++);
    _logging->debug("W = %d\n", _W);
}

void LucasUVMulFast::release()
{
    _recovery_op = 0;
    _X.reset();
    _R.reset();
    _D.reset();
    _arithmetic.reset();
}

void LucasUVMulFast::execute()
{
    if (result() != nullptr)
        return;

    int i, j;
    int pos, first;
    Giant tmp;
    Giant s;
    std::vector<int16_t> naf_w;
    std::unique_ptr<LucasUV> X1;
//    Giant iX, iD, iR;
//#define DEBUG_INDEX(x) x
#define DEBUG_INDEX(x)

    arithmetic().set_gw(gw().carefully());
    LucasVArithmetic lucas(gw(), _negativeQ);
    LucasV Vn(lucas, 2, false);
    LucasV Vn1(lucas, _P, true);
    DEBUG_INDEX(iX = 0);

    first = iterations()%_L;
    if (first == 0)
        first = _L;
    i = 0;
    if (state() == nullptr)
    {
        DEBUG_INDEX(iR = 0);
        arithmetic().init(R());
        State* tmp_state = new State();
        tmp_state->set(0, Vn, Vn1);
        init_state(tmp_state);
        state()->set_written();
    }
    else
    {
        i = state()->iteration();
        state()->to_Lucas(Vn, Vn1);
        DEBUG_INDEX(iR = iX);
        arithmetic().init(Vn, Vn1, R());
    }
    if (state_check() == nullptr)
    {
        DEBUG_INDEX(iD = 0);
        arithmetic().init(D());
    }
    else
    {
        i = state_check()->iteration();
        state_check()->to_Lucas(Vn, Vn1, D());
    }
    pos = iterations() - i;
    pos = pos + (_L - pos%_L)%_L;
    s = 0;
    for (j = pos + (_L2 - pos%_L2)%_L2; j > pos; j -= _L)
    {
        _exp.arithmetic().substr(_exp, j - _L, _L, tmp);
        s += tmp;
    }
    bool Dinit = (i >= first && pos%_L2 != 0);

    lucas.set_gw(i < 30 || pos <= 30 ? gw().carefully() : gw());

    while (i < iterations())
    {
        pos -= _L;
        _exp.arithmetic().substr(_exp, pos, _L, tmp);
        s += tmp;

        for (j = (i < first ? first - i : _L - (i - first)%_L) - 1; j >= 0; j--, i++, commit_execute<StrongCheckState>(i, state()->iteration(), Vn, Vn1, D()))
        {
            if (i == 30)
                lucas.set_gw(gw());
            if (pos + j == 30)
                lucas.set_gw(gw().carefully());
            if (tmp.bit(j))
            {
                DEBUG_INDEX(iX <<= 1);
                DEBUG_INDEX(iX += 1);
                lucas.add(Vn1, Vn, _P, Vn, GWMUL_STARTNEXTFFT_IF(!is_last(i)));
                lucas.dbl(Vn1, Vn1, GWMUL_STARTNEXTFFT_IF(!is_last(i) && j > 0));
            }
            else
            {
                lucas.add(Vn1, Vn, _P, Vn1, GWMUL_STARTNEXTFFT_IF(!is_last(i) && j > 0));
                DEBUG_INDEX(iX <<= 1);
                lucas.dbl(Vn, Vn, GWMUL_STARTNEXTFFT_IF(!is_last(i)));
            }
            if (j > 0)
                continue;

            arithmetic().set_gw(lucas.gw());
            arithmetic().init(Vn, Vn1, X());
            if (pos%_L2 == 0)
                break;
            if (Dinit)
            {
                bool last = ((pos - _L)%_L2 == 0);
                DEBUG_INDEX(iD += iX);
                arithmetic().add(D(), X(), D(), GWMUL_STARTNEXTFFT_IF(!is_last(i) && !last) | (last ? LucasUVArithmetic::LUCASADD_OPTIMIZE : 0));
            }
            else
            {
                DEBUG_INDEX(iD = iX);
                D() = X();
            }
            Dinit = true;
        }

        if (pos%_L2 == 0)
        {
            i++;
            check();
            _logging->progress().update(i/(double)iterations(), ops());
            if (!_tmp_state_recovery)
                _tmp_state_recovery.reset(new State());
            _tmp_state_recovery->set(i, Vn, Vn1);
            DEBUG_INDEX(tmp = iX);

            _logging->debug("performing Gerbicz-Li check at %d,%d.\n", i, pos);
            arithmetic().set_gw(gw().carefully());

            naf_w.clear();
            get_NAF_W(_W, s, naf_w, false);
            s = 0;
            if (Dinit)
            {
                DEBUG_INDEX(iR += iD);
                arithmetic().add(R(), D(), R());
            }
            if (naf_w.size() <= _L)
            {
                DEBUG_INDEX(swap(iX, iR));
                swap(X(), R());
            }
            for (j = (naf_w.size() > _L ? (int)naf_w.size() : _L) - 1; j >= 0; j--)
            {
                if (j >= _L && j == naf_w.size() - 1)
                {
                    DEBUG_INDEX(iX = naf_w[j]);
                    arithmetic().init_small(naf_w[j], X());
                }
                else if (j < naf_w.size() && naf_w[j] != 0)
                {
                    DEBUG_INDEX(iX <<= 1);
                    DEBUG_INDEX(iX += naf_w[j]);
                    arithmetic().dbl_add_small(X(), naf_w[j], X(), 0);
                }
                else
                {
                    DEBUG_INDEX(iX <<= 1);
                    arithmetic().dbl(X(), X(), 0);
                }
                if (j == _L)
                {
                    DEBUG_INDEX(iX += iR);
                    arithmetic().add(X(), R(), X());
                }
            }
            DEBUG_INDEX(swap(iX, iR));
            swap(X(), R());

            DEBUG_INDEX(iX = tmp);
            _tmp_state_recovery->to_Lucas(Vn, Vn1);
            arithmetic().init(Vn, Vn1, X());
            if (Dinit)
            {
                DEBUG_INDEX(iD += iX);
                arithmetic().add(D(), X(), D());
            }
            else
            {
                DEBUG_INDEX(iD = iX);
                D() = X();
            }
            DEBUG_INDEX(iD -= iR);
            DEBUG_INDEX(GWASSERT(iD == 0));
            gw().carefully().sub(D().V(), R().V(), D().V(), 0);
            if (D().V() != 0 || (R().U() == 0 && R().V() == 0))
            {
                _logging->error("Gerbicz-Li check failed at %.1f%%.\n", 100.0*i/iterations());
                if (_file != nullptr)
                    _file->clear();
                _state.reset(new TaskState(5));
                _state->set(_state_recovery->iteration());
                _restart_op = _recovery_op;
                throw TaskRestartException();
            }

            DEBUG_INDEX(iR = iX);
            swap(X(), R());
            arithmetic().init(D());
            Dinit = false;
            _tmp_state_recovery.swap(_state_recovery);
            _state.reset(new TaskState(5));
            _state->set(_state_recovery->iteration());
            _logging->progress().update(_state_recovery->iteration()/(double)iterations(), ops());
            on_state();
            _recovery_op = _restart_op;
            _restart_count = 0;
        }
    }

    done();
}
