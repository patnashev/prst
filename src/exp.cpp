
#include <cmath>
#include <iostream>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "exp.h"
#include "exception.h"

using namespace arithmetic;

void BaseExp::done()
{
    InputTask::done();
    _logging->set_prefix("");
}

void CarefulExp::init(InputNum* input, GWState* gwstate, Logging* logging)
{
    GWASSERT(!smooth());
    GWASSERT(_x0 != 0 || !_X0.empty());
    GWASSERT(_x0 <= (uint32_t)gwstate->maxmulbyconst);
    BaseExp::init(input, gwstate, nullptr, nullptr, logging, _exp.bitlen() - 1 + (!_tail.empty() ? 1 : 0));
    _state_update_period = MULS_PER_STATE_UPDATE*2/3;
    _logging->set_prefix(input->display_text() + " ");
    if (state() != nullptr)
        _logging->info("restarting at %.1f%%.\n", 100.0*state()->iteration()/iterations());
    if (_exp == 1 && _x0 > 0)
        set_state<State>(0, _x0);
    if (_exp == 1 && !_X0.empty())
        set_state<State>(0, _X0);
}

void CarefulExp::execute()
{
    int i, len;

    len = _exp.bitlen() - 1;
    GWNum X0(gw());
    if (!_X0.empty())
        X0 = _X0;
    if (_x0 > 0)
        gw().setmulbyconst(_x0);
    GWNum X(gw());
    if (state() == nullptr)
    {
        i = 0;
        if (!_X0.empty())
            X = X0;
        if (_x0 > 0)
            X = _x0;
    }
    else
    {
        i = state()->iteration();
        X = state()->X();
    }
    for (; i < len; i++, commit_execute<State>(i, X))
    {
        gw().carefully().square(X, X, (_x0 > 0 && _exp.bit(len - i - 1) ? GWMUL_MULBYCONST : 0));
        if (!_X0.empty() && _exp.bit(len - i - 1))
            gw().carefully().mul(X, X0, X, 0);
    }
    if (i < iterations())
    {
        GWNum T(gw());
        T = _tail;
        gw().carefully().mul(T, X, X, 0);
        i++;
        commit_execute<State>(i, X);
    }

    done();
}

void MultipointExp::init(InputNum* input, GWState* gwstate, File* file, Logging* logging)
{
    GWASSERT(smooth() || _x0 != 0 || !_X0.empty());
    GWASSERT(!smooth() || (_x0 == 0 && _X0.empty()));
    GWASSERT(_x0 <= (uint32_t)gwstate->maxmulbyconst);
    BaseExp::init(input, gwstate, file, nullptr, logging, _points.back() + (!_tail.empty() ? 1 : 0));
    _state_update_period = MULS_PER_STATE_UPDATE;
    if (smooth() && b() != 2)
        _state_update_period = (int)(_state_update_period/log2(b()));
    if (_error_check)
        _logging->info("max roundoff check enabled.\n");
    State* state = read_state<State>(file);
    if (state != nullptr)
        init_state(state);
    if (_points.back() == 0 && _x0 > 0)
        set_state<State>(0, _x0);
    if (_points.back() == 0 && !_X0.empty())
        set_state<State>(0, _X0);
}

void MultipointExp::init_state(State* state)
{
    if (state == nullptr)
    {
        _state.reset();
        return;
    }
    _state.reset(state);
    _logging->progress().update(0, (int)_gwstate->handle.fft_count/2);
    _logging->set_prefix(_input->display_text() + " ");
    if (_state->iteration() > 0)
        _logging->info("restarting at %.1f%%.\n", 100.0*_state->iteration()/iterations());
}

void MultipointExp::setup()
{
    if (!_X)
        _X.reset(new GWNum(gw()));
    if (!smooth() && !_X0.empty() && _U.empty())
    {
        X() = _X0;
        slide_init(_exp.bitlen() - 1);
        std::vector<arithmetic::GWNum> U;
        U.swap(_U);
        commit_setup();
        U.swap(_U);
    }
}

void MultipointExp::release()
{
    _X.reset();
    _U.clear();
}

void MultipointExp::execute()
{
    int i, next_point;
    int len;
    Giant tmp;
    int last_power = -1;

    len = _exp.bitlen() - 1;
    if (_x0 > 0)
        gw().setmulbyconst(_x0);
    if (state() == nullptr && !smooth())
    {
        i = 0;
        if (!_X0.empty())
            X() = _X0;
        if (_x0 > 0)
            X() = _x0;
        tmp = _x0;
        init_state(new State(0, !_X0.empty() ? _X0 : tmp));
        state()->set_written();
    }
    else
    {
        GWASSERT(state() != nullptr);
        i = state()->iteration();
        X() = state()->X();
    }
    if (i < 30)
        gwset_carefully_count(gw().gwdata(), 30 - i);

    for (next_point = 0; next_point < _points.size() && i >= _points[next_point]; next_point++);
    for (; next_point < _points.size(); next_point++)
    {
        if ((smooth() && b() == 2) || (!smooth() && _x0 > 0))
        {
            for (; i < _points[next_point]; i++)
            {
                gw().square(X(), X(), (!smooth() && _exp.bit(len - i - 1) ? GWMUL_MULBYCONST : 0) | GWMUL_STARTNEXTFFT_IF(!is_last(i) && i + 1 != _points[next_point]));
                if (i + 1 != _points[next_point])
                    commit_execute<State>(i + 1, X());
            }
        }
        else if (smooth())
        {
            if (last_power != _points[next_point] - i)
            {
                last_power = _points[next_point] - i;
                tmp = b();
                tmp.power(last_power);
            }
            sliding_window(tmp);
            i = _points[next_point];
        }
        else if (!_X0.empty())
        {
            slide(_exp, i, _points[next_point], true);
            i = _points[next_point];
        }

        check();
        if (!_tmp_state)
            _tmp_state.reset(new State());
        static_cast<State*>(_tmp_state.get())->set(i, X());
        if (_on_point != nullptr)
        {
            _logging->progress().update(_points[next_point]/(double)iterations(), (int)_gwstate->handle.fft_count/2);
            _logging->progress_save();
            if (_on_point(next_point, static_cast<State*>(_tmp_state.get())->X()))
            {
                _tmp_state->set_written();
                _last_write = std::chrono::system_clock::now();
            }
        }
        _tmp_state.swap(_state);
        on_state();
    }
    if (i < iterations())
    {
        GWNum T(gw());
        T = _tail;
        gw().carefully().mul(T, X(), X(), 0);
        i++;
        commit_execute<State>(i, X());
    }

    done();
}

int MultipointExp::slide_init(int len)
{
    int W;
    for (W = 1; (W < _W || _W == -1) && ((1 << (W + 1)) <= _max_size || _max_size == -1) && (1 << (W - 1)) + len*(1 + 1/(W + 1.0)) > (1 << (W - 0)) + len*(1 + 1/(W + 2.0)); W++);

    _U.reserve((size_t)1 << (W - 1));
    if (_U.size() <= 0)
        _U.emplace_back(gw());
    swap(_U[0], X());
    if (W > 1)
        gw().square(_U[0], X(), GWMUL_STARTNEXTFFT);
    for (int i = 1; i < (1 << (W - 1)); i++)
    {
        if (_U.size() <= i)
            _U.emplace_back(gw());
        gw().mul(X(), _U[i - 1], _U[i], GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
    }
    return W;
}

void MultipointExp::slide(const arithmetic::Giant& exp, int start, int end, bool commit, int W)
{
    int i, j;
    int len = exp.bitlen() - 1;
    if (W == 0)
        for (W = 1; ((size_t)1 << W) <= _U.size(); W++);

    for (i = len - start - 1; i >= len - end;)
    {
        if (exp.bit(i) == 0)
        {
            gw().square(X(), X(), GWMUL_STARTNEXTFFT_IF(i > len - end));
            i--;
        }
        else
        {
            j = i - W + 1;
            if (j < len - end)
                j = len - end;
            for (; exp.bit(j) == 0; j++);
            int ui = 0;

            while (i >= j)
            {
                gw().square(X(), X(), GWMUL_STARTNEXTFFT);
                ui <<= 1;
                ui += exp.bit(i) ? 1 : 0;
                i--;
            }

            gw().mul(_U[ui/2], X(), X(), GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(i > len - end));
        }
        if (commit && i >= len - end)
            commit_execute<State>(len - i - 1, X());
    }
}

void MultipointExp::sliding_window(const arithmetic::Giant& exp)
{
    int len = exp.bitlen() - 1;
    int W = slide_init(len);

    int i = len;
    int j = i - W + 1;
    if (j < 0)
        j = 0;
    for (; exp.bit(j) == 0; j++);
    int ui = 0;
    while (i >= j)
    {
        ui <<= 1;
        ui += exp.bit(i) ? 1 : 0;
        i--;
    }
    X() = _U[ui/2];

    slide(exp, len - i - 1, len, false, W);
}

double MultipointExp::cost()
{
    bool simple = dynamic_cast<SlidingWindowExp*>(this) == nullptr;
    if ((smooth() && b() == 2) || (!smooth() && simple))
        return _points.back();
    else if (smooth())
    {
        double log2b = log2(b());
        int W;
        int first = 0;
        if (_points[0] == 0)
            first = 1;
        for (W = 2; (W < _W || _W == -1) && ((1 << (W + 1)) <= _max_size || _max_size == -1) && (1 << (W - 1)) + log2b*_points[first]*(1 + 1/(W + 1.0)) >(1 << (W - 0)) + log2b*_points[first]*(1 + 1/(W + 2.0)); W++);
        double cost = ((1 << (W - 1)) + log2b*_points[first]*(1 + 1/(W + 1.0)));
        if (_points.size() > 1 + first)
        {
            cost *= (_points.size() - 1 - first);
            int last = _points[_points.size() - 1] - _points[_points.size() - 2];
            for (W = 2; (W < _W || _W == -1) && ((1 << (W + 1)) <= _max_size || _max_size == -1) && (1 << (W - 1)) + log2b*last*(1 + 1/(W + 1.0)) >(1 << (W - 0)) + log2b*last*(1 + 1/(W + 2.0)); W++);
            cost += ((1 << (W - 1)) + log2b*last*(1 + 1/(W + 1.0)));
        }
        return cost;
    }
    else
    {
        int len = _exp.bitlen() - 1;
        int W;
        for (W = 2; (W < _W || _W == -1) && ((1 << (W + 1)) <= _max_size || _max_size == -1) && (1 << (W - 1)) + len*(1 + 1/(W + 1.0)) >(1 << (W - 0)) + len*(1 + 1/(W + 2.0)); W++);
        return (1 << (W - 1)) + _points.back()*(1 + 1/(W + 1.0));
    }
}

void StrongCheckMultipointExp::Gerbicz_params(int iters, double log2b, int& L, int &L2)
{
    int i;
    //if (log2b > 1.5)
    //    log2b /= 2;
    log2b = 1;
    L = (int)std::sqrt(iters/log2b);
    L2 = iters - iters%L;
    L = L2/L;
    for (i = L + 1; i*i < 2*iters/log2b; i++)
        if (L2 < iters - iters%i)
        {
            L = i;
            L2 = iters - iters%i;
        }
}

double StrongCheckMultipointExp::cost()
{
    int n = _points.back();
    bool simple = dynamic_cast<LiCheckExp*>(this) == nullptr || dynamic_cast<FastLiCheckExp*>(this) != nullptr;
    if ((smooth() && b() == 2) || (!smooth() && simple))
        return n + n/_L + n/_L2*(_L + (!smooth() ? (_L + std::log2(_L2/_L)) : 0));
    else if (smooth())
    {
        double log2b = log2(b());
        int W;
        for (W = 2; (W < _W || _W == -1) && ((1 << (W + 1)) <= _max_size || _max_size == -1) && (1 << (W - 1)) + log2b*_L*(1 + 1/(W + 1.0)) >(1 << (W - 0)) + log2b*_L*(1 + 1/(W + 2.0)); W++);
        return n/_L + (n/_L + n/_L2)*((1 << (W - 1)) + log2b*_L*(1 + 1/(W + 1.0)));
    }
    else
    {
        int len = _exp.bitlen() - 1;
        int W;
        for (W = 2; (W < _W || _W == -1) && ((1 << (W + 1)) <= _max_size || _max_size == -1) && (1 << (W - 1)) + len*(1 + 1/(W + 1.0)) >(1 << (W - 0)) + len*(1 + 1/(W + 2.0)); W++);
        return (1 << (W - 1)) + n*(1 + 1/(W + 1.0)) + n/_L + n/_L2*(_L + (_L + std::log2(_L2/_L))*(1 + 1/(W + 1.0)));
    }
}

void StrongCheckMultipointExp::init(InputNum* input, GWState* gwstate, File* file, File* file_recovery, Logging* logging)
{
    GWASSERT(smooth() || _x0 != 0 || !_X0.empty());
    GWASSERT(!smooth() || (_x0 == 0 && _X0.empty()));
    GWASSERT(_x0 <= (uint32_t)gwstate->maxmulbyconst);
    GWASSERT(_points.back() > 0);
    BaseExp::init(input, gwstate, file, read_state<StrongCheckState>(file), logging, _points.back() + (!_tail.empty() ? 1 : 0));
    _state_update_period = MULS_PER_STATE_UPDATE;
    if (smooth() && b() != 2)
        _state_update_period = (int)(_state_update_period/log2(b()));
    _logging->info("Gerbicz%s check enabled, L2 = %d*%d.\n", !smooth() ? "-Li" : "", _L, _L2/_L);
    _logging->report_param("L", _L);
    _logging->report_param("L2", _L2);
    if (_error_check)
        _logging->info("max roundoff check enabled.\n");
    _file_recovery = file_recovery;
    _state_recovery.reset();
    State* state_recovery = read_state<State>(file_recovery);
    if (state_recovery != nullptr)
        init_state(state_recovery);
    _file_recovery_empty = state_recovery == nullptr;
}

void StrongCheckMultipointExp::init_state(State* state)
{
    if (state == nullptr)
    {
        _state_recovery.reset();
        _state.reset();
        return;
    }
    _logging->progress().update(0, (int)_gwstate->handle.fft_count/2);
    _logging->set_prefix(_input->display_text() + " ");
    _state_recovery.reset(state);
    _recovery_op = 0;
    if (!_state || (state_check() != nullptr ? state_check()->recovery() : _state->iteration()) != _state_recovery->iteration())
    {
        _state.reset(new TaskState(5));
        _state->set(_state_recovery->iteration());
    }
    if (_state->iteration() > 0)
        _logging->info("restarting at %.1f%%.\n", 100.0*_state->iteration()/iterations());
}

void StrongCheckMultipointExp::StrongCheckState::set(int iteration, int recovery, arithmetic::GWNum& X, arithmetic::GWNum& D)
{
    TaskState::set(iteration);
    _recovery = recovery;
    if (!_gwX)
        _gwX.reset(new GWNum(X.arithmetic()));
    *_gwX = X;
    if (!_gwD)
        _gwD.reset(new GWNum(D.arithmetic()));
    *_gwD = D;
}

void StrongCheckMultipointExp::write_state()
{
    if (_file_recovery != nullptr && _state_recovery && !_state_recovery->is_written())
    {
        _file_recovery->write(*_state_recovery);
        _file_recovery_empty = false;
    }
    if (state_check() != nullptr)
    {
        try
        {
            if (state_check()->gwX())
                state_check()->X() = *state_check()->gwX();
            if (state_check()->gwD())
                state_check()->D() = *state_check()->gwD();
        }
        catch (const ArithmeticException&)
        {
            _state.reset(new TaskState(5));
            _state->set(_state_recovery->iteration());
            throw;
        }
    }
    Task::write_state();
}

void StrongCheckMultipointExp::setup()
{
    MultipointExp::setup();
    if (!_R)
        _R.reset(new GWNum(gw()));
    if (!_D)
        _D.reset(new GWNum(gw()));
    if (state_check() != nullptr)
    {
        if (!state_check()->gwX())
        {
            state_check()->gwX().reset(new GWNum(gw()));
            *state_check()->gwX() = state_check()->X();
        }
        if (!state_check()->gwD())
        {
            state_check()->gwD().reset(new GWNum(gw()));
            *state_check()->gwD() = state_check()->D();
        }
    }
}

void StrongCheckMultipointExp::release()
{
    _recovery_op = 0;
    _R.reset();
    _D.reset();
    if (state_check() != nullptr)
    {
        try
        {
            if (state_check()->gwX())
                state_check()->X() = *state_check()->gwX();
            state_check()->gwX().reset();
            if (state_check()->gwD())
                state_check()->D() = *state_check()->gwD();
            state_check()->gwD().reset();
        }
        catch (const ArithmeticException&)
        {
            _state.reset(new TaskState(5));
            _state->set(_state_recovery->iteration());
        }
    }
    _tmp_state.reset();
    MultipointExp::release();
}

void StrongCheckMultipointExp::execute()
{
    int i, j, next_point, next_check;
    int len, l;
    Giant exp;
    int last_power = -1;
    Giant tmp;
    Giant tmp2;

    len = _exp.bitlen() - 1;
    GWNum X0(gw());
    if (!_X0.empty())
        X0 = _X0;
    if (_x0 > 0)
        gw().setmulbyconst(_x0);
    if (state() == nullptr && !smooth())
    {
        i = 0;
        if (!_X0.empty())
            R() = X0;
        if (_x0 > 0)
            R() = _x0;
        tmp = _x0;
        init_state(new State(0, !_X0.empty() ? _X0 : tmp));
        state()->set_written();
    }
    else
    {
        GWASSERT(state() != nullptr);
        i = state()->iteration();
        R() = state()->X();
    }
    if (state_check() == nullptr)
    {
        X() = R();
        D() = R();
    }
    else
    {
        i = state_check()->iteration();
        X() = *state_check()->gwX();
        D() = *state_check()->gwD();
    }
    if (i < 30)
        gwset_carefully_count(gw().gwdata(), 30 - i);

    for (next_point = 0; next_point < _points.size() && state()->iteration() >= abs(_points[next_point]); next_point++);
    while (next_point < _points.size())
    {
        for (next_check = next_point; _points[next_check] < 0; next_check++);
        int L = _L;
        int L2 = _L2;
        while ((_points[next_check] - state()->iteration()) < L2 && L > 1)
        {
            if (L == 3)
                L = 2;
            else
                L /= 2;
            L2 = L*L;
            last_power = -1;
            for (next_check = next_point; _points[next_check] < 0; next_check++);
        }
        if (i - state()->iteration() > L2)
        {
            i = state()->iteration();
            X() = R();
            D() = R();
            _state.reset(new TaskState(5));
            _state->set(i);
        }
        else
            for (; next_point < next_check && i >= abs(_points[next_point]); next_point++);

        if ((smooth() && b() == 2) || (!smooth() && _x0 > 0))
        {
            for (j = i - state()->iteration(); j < L2; j++, i++, commit_execute<StrongCheckState>(i, state()->iteration(), X(), D()))
            {
                gw().square(X(), X(), (!smooth() && _exp.bit(len - i - 1) ? GWMUL_MULBYCONST : 0) | GWMUL_STARTNEXTFFT_IF(!is_last(i) && i + 1 != abs(_points[next_point]) && j + 1 != L2));
                if (j + 1 != L2 && i + 1 == -_points[next_point])
                {
                    check();
                    _logging->progress().update(-_points[next_point]/(double)iterations(), (int)_gwstate->handle.fft_count/2);
                    _logging->progress_save();
                    tmp = X();
                    if (_on_point != nullptr)
                        _on_point(next_point, tmp);
                    _logging->progress().update(-_points[next_point]/(double)iterations(), (int)_gwstate->handle.fft_count/2);
                    set_state<StrongCheckState>(i + 1, state()->iteration(), X(), D());
                    next_point++;
                }
                if (j + 1 != L2 && (j + 1)%L == 0)
                    gw().mul(X(), D(), D(), is_last(i) ? GWMUL_PRESERVE_S1 : GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(j + 1 + L != L2));
            }
        }
        else if (smooth())
        {
            GWASSERT((i - state()->iteration())%L == 0);
            for (j = i - state()->iteration(); j < L2; j += L, i += L, commit_execute<StrongCheckState>(i, state()->iteration(), X(), D()))
            {
                if (last_power != L)
                {
                    last_power = L;
                    exp = b();
                    exp.power(last_power);
                }
                sliding_window(exp);
                if (j + L != L2 && i + L == -_points[next_point])
                {
                    check();
                    _logging->progress().update(-_points[next_point]/(double)iterations(), (int)_gwstate->handle.fft_count/2);
                    _logging->progress_save();
                    tmp = X();
                    if (_on_point != nullptr)
                        _on_point(next_point, tmp);
                    _logging->progress().update(-_points[next_point]/(double)iterations(), (int)_gwstate->handle.fft_count/2);
                    set_state<StrongCheckState>(i + L, state()->iteration(), X(), D());
                    next_point++;
                }
                if (j + L != L2)
                    gw().mul(X(), D(), D(), is_last(i + L - 1) ? GWMUL_PRESERVE_S1 : GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(j + L + L != L2));
            }
        }
        else if(!_X0.empty())
        {
            GWASSERT((i - state()->iteration())%L == 0);
            for (j = i - state()->iteration(); j < L2; j += L, i += L, commit_execute<StrongCheckState>(i, state()->iteration(), X(), D()))
            {
                if (_points[next_point] < 0 && i + L >= -_points[next_point])
                {
                    slide(_exp, i, -_points[next_point], false);
                    check();
                    _logging->progress().update(-_points[next_point]/(double)iterations(), (int)_gwstate->handle.fft_count/2);
                    _logging->progress_save();
                    tmp = X();
                    if (_on_point != nullptr)
                        _on_point(next_point, tmp);
                    _logging->progress().update(-_points[next_point]/(double)iterations(), (int)_gwstate->handle.fft_count/2);
                    set_state<StrongCheckState>(-_points[next_point], state()->iteration(), X(), D());
                    slide(_exp, -_points[next_point], i + L, false);
                    next_point++;
                }
                else
                    slide(_exp, i, i + L, false);
                if (j + L != L2)
                    gw().mul(X(), D(), D(), is_last(i + L - 1) ? GWMUL_PRESERVE_S1 : GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(j + L + L != L2));
            }
        }
        check();
        _logging->progress().update(i/(double)iterations(), (int)_gwstate->handle.fft_count/2);
        if (next_point != next_check)
        {
            _logging->error("point missed, invalid parameters.\n");
            throw TaskAbortException();
        }

        _logging->debug("performing Gerbicz%s check at %d,%d, L2 = %d*%d.\n", !smooth() ? "-Li" : "", next_check, i, L, L2/L);
        GWNum T(D());
        gw().carefully().mul(X(), D(), D(), 0);
        swap(T, X());
        if ((smooth() && b() == 2) || !smooth())
        {
            for (j = 0; j < L; j++)
                gw().carefully().square(X(), X(), 0);
        }
        else if (smooth())
        {
            GWArithmetic* tmpgw = _gw;
            _gw = &gw().carefully();
            if (last_power != L)
            {
                last_power = L;
                exp = b();
                exp.power(last_power);
            }
            sliding_window(exp);
            _gw = tmpgw;
        }
        if (!smooth())
        {
            for (j = 0, tmp = 0; j*L < L2; j++, tmp += tmp2)
                _exp.arithmetic().substr(_exp, len - state()->iteration() - j*L - L, L, tmp2);
            if (tmp != 0)
            {
                GWNum TX(gw());
                if (!_X0.empty())
                {
                    TX = X0;
                    swap(TX, X());
                    //GWArithmetic* tmpgw = _gw;
                    //_gw = &gw().carefully();
                    slide(tmp, 0, tmp.bitlen() - 1, false);
                    //_gw = tmpgw;
                    swap(TX, X());
                }
                if (_x0 > 0)
                {
                    TX = _x0;
                    for (l = tmp.bitlen() - 1, j = 0; j < l; j++)
                        gw().carefully().square(TX, TX, tmp.bit(l - j - 1) ? GWMUL_MULBYCONST : 0);
                }
                gw().carefully().mul(TX, X(), X(), 0);
            }
        }
        gw().carefully().mul(R(), X(), X(), 0);
        gw().carefully().sub(X(), D(), X(), 0);
        swap(T, X());
        if (T != 0 || D() == 0)
        {
            _logging->error("Gerbicz%s check failed at %.1f%%.\n", !smooth() ? "-Li" : "", 100.0*i/iterations());
            if (_file != nullptr)
                _file->clear();
            _state.reset(new TaskState(5));
            _state->set(_state_recovery->iteration());
            _restart_op = _recovery_op;
            throw TaskRestartException();
        }

        R() = X();
        D() = X();
        if (!_tmp_state_recovery)
            _tmp_state_recovery.reset(new State());
        _tmp_state_recovery->set(i, R());
        if (i == _points[next_check])
        {
            if (_on_point != nullptr)
            {
                _logging->progress_save();
                if (_on_point(next_check, _tmp_state_recovery->X()))
                {
                    _tmp_state_recovery->set_written();
                    _last_write = std::chrono::system_clock::now();
                    if (!_file_recovery_empty)
                    {
                        _file_recovery->clear();
                        _file_recovery_empty = true;
                    }
                }
            }
            next_point++;
        }
        _tmp_state_recovery.swap(_state_recovery);
        _state.reset(new TaskState(5));
        _state->set(i);
        _logging->progress().update(i/(double)iterations(), (int)_gwstate->handle.fft_count/2);
        on_state();
        _recovery_op = _restart_op;
        _restart_count = 0;
    }
    if (i < iterations())
    {
        D() = _tail;
        gw().carefully().mul(D(), R(), R(), 0);
        i++;
        tmp = R();
        _state_recovery.reset(new State(i, std::move(tmp)));
        _state->set(i);
        on_state();
    }

    done();
}
