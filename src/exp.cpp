
#include <cmath>
#include <iostream>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "exp.h"
#include "exception.h"

using namespace arithmetic;

BaseExp::State* BaseExp::State::read_file(File* file)
{
    if (file == nullptr)
        return nullptr;
    std::unique_ptr<Reader> reader(file->get_reader());
    if (!reader)
        return nullptr;
    std::unique_ptr<State> state;
    if (reader->type() == StateValue::TYPE)
        state.reset(new StateValue());
    else if (reader->type() == StateSerialized::TYPE)
        state.reset(new StateSerialized());
    else
        return nullptr;
    if (!state->read(*reader))
        return nullptr;
    state->set_written();
    return state.release();
}

BaseExp::State* BaseExp::State::read_file(File* file, StateValue* value, StateSerialized* serialized)
{
    if (file == nullptr)
        return nullptr;
    std::unique_ptr<Reader> reader(file->get_reader());
    if (!reader)
        return nullptr;
    State* state;
    if (reader->type() == StateValue::TYPE)
        state = value;
    else if (reader->type() == StateSerialized::TYPE)
        state = serialized;
    else
        return nullptr;
    if (!state->read(*reader))
        return nullptr;
    state->set_written();
    return state;
}

BaseExp::State* BaseExp::State::cast(bool value, std::unique_ptr<TaskState>& state)
{
    State* cast_state;
    if (value)
    {
        if (!state || (cast_state = dynamic_cast<StateValue*>(state.get())) == nullptr)
            state.reset(cast_state = new StateValue());
    }
    else
    {
        if (!state || (cast_state = dynamic_cast<StateSerialized*>(state.get())) == nullptr)
            state.reset(cast_state = new StateSerialized());
    }
    return cast_state;
}

void CarefulExp::init(InputNum* input, GWState* gwstate, Logging* logging)
{
    GWASSERT(!smooth());
    GWASSERT(_x0 != 0 || !_X0.empty());
    GWASSERT(_x0 <= (uint32_t)gwstate->maxmulbyconst);
    BaseExp::init(input, gwstate, nullptr, nullptr, logging, _exp.bitlen() - 1 + (!_tail.empty() ? 1 : 0));
    _state_update_period = MULS_PER_STATE_UPDATE*2/3;
    if (_exp == 1 && _x0 > 0)
        set_state<StateValue>(0, _x0);
    if (_exp == 1 && !_X0.empty())
        set_state<StateValue>(0, _X0);
}

void CarefulExp::execute()
{
    if (result() != nullptr)
        return;

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
        state()->to_GWNum(X);
    }
    for (; i < len; i++, commit_execute(i, X))
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
        commit_execute(i, X);
    }

    done();
}

void MultipointExp::init(InputNum* input, GWState* gwstate, File* file, Logging* logging)
{
    GWASSERT(smooth() || _x0 != 0 || !_X0.empty());
    GWASSERT(!smooth() || (_x0 == 0 && _X0.empty()));
    GWASSERT(_x0 <= (uint32_t)gwstate->maxmulbyconst);
    BaseExp::init(input, gwstate, file, nullptr, logging, _points.back().pos + (!_tail.empty() ? 1 : 0));
    _state_update_period = MULS_PER_STATE_UPDATE;
    if (smooth() && b() != 2)
        _state_update_period = (int)(_state_update_period/log2(b()));
    if (_error_check)
        _logging->info("max roundoff check enabled.\n");
    State* state = State::read_file(file);
    if (state != nullptr)
        init_state(state);
    if (_points.back().pos == 0 && _x0 > 0)
        set_state<StateValue>(0, _x0);
    if (_points.back().pos == 0 && !_X0.empty())
        set_state<StateValue>(0, _X0);
}

void MultipointExp::init_state(State* state)
{
    if (state == nullptr)
    {
        _state.reset();
        return;
    }
    _state.reset(state);
    _logging->progress().update(progress(), 0);
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
        gwset_carefully_count(gw().gwdata(), 30);
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
    if (result() != nullptr)
        return;

    int i, next_point;
    int len;
    Giant tmp;
    State* tmp_state;
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
        tmp_state = new StateSerialized();
        tmp_state->set(0, X());
        init_state(tmp_state);
        state()->set_written();
    }
    else
    {
        GWASSERT(state() != nullptr);
        i = state()->iteration();
        state()->to_GWNum(X());
    }
    if (i < 30)
        gwset_carefully_count(gw().gwdata(), 30 - i);

    for (next_point = 0; next_point < _points.size() && i >= _points[next_point].pos; next_point++);
    for (; next_point < _points.size(); next_point++)
    {
        if ((smooth() && b() == 2) || (!smooth() && _x0 > 0))
        {
            for (; i < _points[next_point].pos; i++)
            {
                (i < iterations() - 30 ? gw() : gw().carefully()).square(X(), X(), (!smooth() && _exp.bit(len - i - 1) ? GWMUL_MULBYCONST : 0) | GWMUL_STARTNEXTFFT_IF(!is_last(i) && i + 1 != _points[next_point].pos));
                if (i + 1 != _points[next_point].pos)
                    commit_execute(i + 1, X());
            }
        }
        else if (smooth())
        {
            if (last_power != _points[next_point].pos - i)
            {
                last_power = _points[next_point].pos - i;
                tmp = b();
                tmp.power(last_power);
            }
            sliding_window(tmp);
            i = _points[next_point].pos;
        }
        else if (!_X0.empty())
        {
            slide(_exp, i, _points[next_point].pos, true);
            i = _points[next_point].pos;
        }

        check();
        tmp_state = State::cast(_points[next_point].value, _tmp_state);
        tmp_state->set(i, X());
        if (_on_point != nullptr)
        {
            _logging->progress().update(_points[next_point].pos/(double)iterations(), ops());
            _logging->progress_save();
            if (_on_point(next_point, tmp_state, *_logging))
            {
                tmp_state->set_written();
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
        commit_execute(i, X());
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
            commit_execute(len - i - 1, X());
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
        return _points.back().pos;
    else if (smooth())
    {
        double log2b = log2(b());
        int W;
        int first = 0;
        if (_points[0].pos == 0)
            first = 1;
        for (W = 2; (W < _W || _W == -1) && ((1 << (W + 1)) <= _max_size || _max_size == -1) && (1 << (W - 1)) + log2b*_points[first].pos*(1 + 1/(W + 1.0)) >(1 << (W - 0)) + log2b*_points[first].pos*(1 + 1/(W + 2.0)); W++);
        double cost = ((1 << (W - 1)) + log2b*_points[first].pos*(1 + 1/(W + 1.0)));
        if (_points.size() > 1 + first)
        {
            cost *= (_points.size() - 1 - first);
            int last = _points[_points.size() - 1].pos - _points[_points.size() - 2].pos;
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
        return (1 << (W - 1)) + _points.back().pos*(1 + 1/(W + 1.0));
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
    int n = _points.back().pos;
    bool simple = dynamic_cast<LiCheckExp*>(this) == nullptr || dynamic_cast<FastLiCheckExp*>(this) != nullptr;
    if ((smooth() && b() == 2) || (!smooth() && simple))
        return n + n/_L + n/_L2*(!smooth() ? _L + std::log2(_L2/_L) : _L)*1.5;
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
        return (1 << (W - 1)) + n*(1 + 1/(W + 1.0)) + n/_L + n/_L2*(_L + std::log2(_L2/_L))*(1 + 1/(W + 1.0));
    }
}

void StrongCheckMultipointExp::init(InputNum* input, GWState* gwstate, File* file, File* file_recovery, Logging* logging)
{
    GWASSERT(smooth() || _x0 != 0 || !_X0.empty());
    GWASSERT(!smooth() || (_x0 == 0 && _X0.empty()));
    GWASSERT(_x0 <= (uint32_t)gwstate->maxmulbyconst);
    GWASSERT(_points.back().pos > 0);
    BaseExp::init(input, gwstate, file, read_state<StrongCheckState>(file), logging, _points.back().pos + (!_tail.empty() ? 1 : 0));
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
    State* state_recovery = State::read_file(file_recovery);
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
}

void StrongCheckMultipointExp::write_state()
{
    if (_file_recovery != nullptr && _state_recovery && !_state_recovery->is_written())
    {
        _file_recovery->write(*_state_recovery);
        _file_recovery_empty = false;
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
}

void StrongCheckMultipointExp::release()
{
    _recovery_op = 0;
    _R.reset();
    _D.reset();
    MultipointExp::release();
}

void StrongCheckMultipointExp::execute()
{
    if (result() != nullptr)
        return;

    int i, j, next_point, next_check;
    int len;
    Giant exp;
    int last_power = -1;
    Giant tmp;
    Giant tmp2;
    State* tmp_state;

    len = _exp.bitlen() - 1;
    GWNum X0(gw());
    if (!_X0.empty())
        X0 = _X0;
    if (_x0 > 0)
        gw().setmulbyconst(_x0);
    if (state() == nullptr && !smooth())
    {
        i = 0;
        StateValue* state0 = nullptr;
        if (!_X0.empty())
            state0 = new StateValue(0, _X0);
        if (_x0 > 0)
            state0 = new StateValue(0, _x0);
        state0->to_GWNum(R());
        tmp = R();
        if (tmp != state0->value())
        {
            _logging->warning("Value initialization error.\n");
            throw TaskRestartException();
        }
        init_state(state0);
        state()->set_written();
    }
    else
    {
        GWASSERT(state() != nullptr);
        i = state()->iteration();
        state()->to_GWNum(R());
    }
    if (state_check() == nullptr)
    {
        X() = R();
        D() = R();
    }
    else
    {
        i = state_check()->iteration();
        X() = state_check()->X();
        D() = state_check()->D();
    }
    if (i < 30)
        gwset_carefully_count(gw().gwdata(), 30 - i);

    for (next_point = 0; next_point < _points.size() && state()->iteration() >= _points[next_point].pos; next_point++);
    while (next_point < _points.size())
    {
        for (next_check = next_point; !_points[next_check].check; next_check++);
        int L = _L;
        int L2 = _L2;
        while ((_points[next_check].pos - state()->iteration()) < L2 && L > 1)
        {
            if (L == 3)
                L = 2;
            else
                L /= 2;
            L2 = L*L;
            last_power = -1;
            for (next_check = next_point; !_points[next_check].check; next_check++);
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
            for (; next_point < next_check && i >= _points[next_point].pos; next_point++);

        if ((smooth() && b() == 2) || (!smooth() && _x0 > 0))
        {
            for (j = i - state()->iteration(); j < L2; j++, i++, Task::commit_execute<StrongCheckState>(i, state()->iteration(), X(), D()))
            {
                (i < iterations() - 30 ? gw() : gw().carefully()).square(X(), X(), (!smooth() && _exp.bit(len - i - 1) ? GWMUL_MULBYCONST : 0) | GWMUL_STARTNEXTFFT_IF(!is_last(i) && i + 1 != _points[next_point].pos && j + 1 != L2));
                if (j + 1 != L2 && i + 1 == _points[next_point].pos && !_points[next_point].check)
                {
                    check();
                    _logging->progress().update(_points[next_point].pos/(double)iterations(), ops());
                    _logging->progress_save();
                    if (_on_point != nullptr)
                    {
                        tmp_state = State::cast(_points[next_point].value, _tmp_state_recovery);
                        tmp_state->set(_points[next_point].pos, X());
                        _on_point(next_point, tmp_state, *_logging);
                    }
                    _logging->progress().update(_points[next_point].pos/(double)iterations(), ops());
                    set_state<StrongCheckState>(_points[next_point].pos, state()->iteration(), X(), D());
                    next_point++;
                }
                if (j + 1 != L2 && (j + 1)%L == 0)
                    gw().mul(X(), D(), D(), is_last(i) ? GWMUL_PRESERVE_S1 : GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(j + 1 + L != L2));
            }
        }
        else if (smooth())
        {
            GWASSERT((i - state()->iteration())%L == 0);
            for (j = i - state()->iteration(); j < L2; j += L, i += L, Task::commit_execute<StrongCheckState>(i, state()->iteration(), X(), D()))
            {
                if (last_power != L)
                {
                    last_power = L;
                    exp = b();
                    exp.power(last_power);
                }
                sliding_window(exp);
                if (j + L != L2 && i + L == _points[next_point].pos && !_points[next_point].check)
                {
                    check();
                    _logging->progress().update(_points[next_point].pos/(double)iterations(), ops());
                    _logging->progress_save();
                    if (_on_point != nullptr)
                    {
                        tmp_state = State::cast(_points[next_point].value, _tmp_state_recovery);
                        tmp_state->set(_points[next_point].pos, X());
                        _on_point(next_point, tmp_state, *_logging);
                    }
                    _logging->progress().update(_points[next_point].pos/(double)iterations(), ops());
                    set_state<StrongCheckState>(_points[next_point].pos, state()->iteration(), X(), D());
                    next_point++;
                }
                if (j + L != L2)
                    gw().mul(X(), D(), D(), is_last(i + L - 1) ? GWMUL_PRESERVE_S1 : GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(j + L + L != L2));
            }
        }
        else if(!_X0.empty())
        {
            GWASSERT((i - state()->iteration())%L == 0);
            for (j = i - state()->iteration(); j < L2; j += L, i += L, Task::commit_execute<StrongCheckState>(i, state()->iteration(), X(), D()))
            {
                if (i + L >= _points[next_point].pos && !_points[next_point].check)
                {
                    slide(_exp, i, _points[next_point].pos, false);
                    check();
                    _logging->progress().update(_points[next_point].pos/(double)iterations(), ops());
                    _logging->progress_save();
                    if (_on_point != nullptr)
                    {
                        tmp_state = State::cast(_points[next_point].value, _tmp_state_recovery);
                        tmp_state->set(_points[next_point].pos, X());
                        _on_point(next_point, tmp_state, *_logging);
                    }
                    _logging->progress().update(_points[next_point].pos/(double)iterations(), ops());
                    set_state<StrongCheckState>(_points[next_point].pos, state()->iteration(), X(), D());
                    slide(_exp, _points[next_point].pos, i + L, false);
                    next_point++;
                }
                else
                    slide(_exp, i, i + L, false);
                if (j + L != L2)
                    gw().mul(X(), D(), D(), is_last(i + L - 1) ? GWMUL_PRESERVE_S1 : GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(j + L + L != L2));
            }
        }
        check();
        _logging->progress().update(i/(double)iterations(), ops());
        if (next_point != next_check)
        {
            _logging->error("point missed, invalid parameters.\n");
            throw TaskAbortException();
        }

        _logging->debug("performing Gerbicz%s check at %d,%d, L2 = %d*%d.\n", !smooth() ? "-Li" : "", next_check, i, L, L2/L);
        tmp_state = State::cast(i == _points[next_check].pos && _points[next_check].value, _tmp_state_recovery);
        tmp_state->set(i, X());
        GWNum T(gw());
        tmp_state->to_GWNum(T);
        gw().carefully().mul(T, D(), X(), 0);
        swap(X(), D());
        if (smooth() && b() == 2)
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
            if (tmp.bitlen() > L)
            {
                tmp.arithmetic().shiftright(tmp, L, tmp2);
                tmp.arithmetic().substr(tmp, 0, L, tmp);
                GWNum TX(gw());
                if (!_X0.empty())
                {
                    TX = X0;
                    swap(TX, X());
                    //GWArithmetic* tmpgw = _gw;
                    //_gw = &gw().carefully();
                    slide(tmp2, 0, tmp2.bitlen() - 1, false);
                    //_gw = tmpgw;
                    swap(TX, X());
                }
                if (_x0 > 0)
                {
                    TX = _x0;
                    for (j = tmp2.bitlen() - 2; j >= 0; j--)
                        gw().carefully().square(TX, TX, tmp2.bit(j) ? GWMUL_MULBYCONST : 0);
                }
                gw().carefully().mul(TX, X(), X(), 0);
            }
            if (!_X0.empty())
            {
                //GWArithmetic* tmpgw = _gw;
                //_gw = &gw().carefully();
                tmp2 = 1;
                tmp2 <<= L;
                tmp += tmp2;
                slide(tmp, 0, L, false);
                //_gw = tmpgw;
            }
            if (_x0 > 0)
            {
                for (j = L - 1; j >= 0; j--)
                    gw().carefully().square(X(), X(), tmp.bit(j) ? GWMUL_MULBYCONST : 0);
            }
        }
        gw().carefully().mul(R(), X(), X(), 0);
        gw().carefully().sub(X(), D(), X(), 0);
        swap(T, X());
        tmp = T;
        if (_gwstate->need_mod() && dynamic_cast<StateValue*>(tmp_state))
            _gwstate->mod(tmp, tmp);
        tmp2 = D();
        if (_gwstate->need_mod() && dynamic_cast<StateValue*>(tmp_state))
            _gwstate->mod(tmp2, tmp2);
        if (tmp != 0 || tmp2 == 0)
        {
            _logging->warning("Gerbicz%s check failed at %.1f%%.\n", !smooth() ? "-Li" : "", 100.0*i/iterations());
            if (_file != nullptr)
                _file->clear();
            _state.reset(new TaskState(5));
            _state->set(_state_recovery->iteration());
            _restart_op = _recovery_op;
            throw TaskRestartException();
        }

        R() = X();
        D() = X();
        if (i == _points[next_check].pos)
        {
            if (_on_point != nullptr)
            {
                _logging->progress_save();
                if (_on_point(next_check, tmp_state, *_logging))
                {
                    tmp_state->set_written();
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
        _logging->progress().update(i/(double)iterations(), ops());
        on_state();
        _recovery_op = _restart_op;
        _restart_count = 0;
    }
    if (i < iterations())
    {
        D() = _tail;
        gw().carefully().mul(D(), R(), R(), 0);
        i++;
        tmp_state = State::cast(true, _tmp_state_recovery);
        tmp_state->set(i, R());
        if (!_tail_inv.empty())
        {
            tmp_state->to_GWNum(R());
            D() = _tail_inv;
            gw().carefully().mul(D(), R(), R(), 0);
            tmp = R();
            if (_gwstate->need_mod())
                _gwstate->mod(tmp, tmp);
            if (tmp != static_cast<StateValue*>(_state_recovery.get())->value())
            {
                R() = _tail;
                gw().carefully().mul(D(), R(), R(), 0);
                tmp2 = R();
                if (_gwstate->need_mod())
                    _gwstate->mod(tmp2, tmp2);
                if (tmp2 != 1)
                {
                    _logging->error("Tail calculation error.\n");
                    throw TaskAbortException();
                }
                else
                {
                    _logging->warning("Tail check failed.\n");
                    throw TaskRestartException();
                }
            }
        }
        _tmp_state_recovery.swap(_state_recovery);
        _state->set(i);
        on_state();
    }

    done();
}
