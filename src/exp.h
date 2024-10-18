#pragma once

#include <functional>
#include "arithmetic.h"
#include "group.h"
#include "integer.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"

class BaseExp : public InputTask
{
public:
    class StateSerialized;
    class StateValue;
    class State : public TaskState
    {
    public:
        State(char type) : TaskState(type) { }
        virtual void set(int iteration, arithmetic::GWNum& X) = 0;
        virtual void to_GWNum(arithmetic::GWNum& X) = 0;
        static State* read_file(File* file);
        static State* read_file(File* file, StateValue* value, StateSerialized* serialized);
        static State* cast(bool value, std::unique_ptr<TaskState>& state);
    };
    class StateSerialized : public State
    {
    public:
        static const char TYPE = 8;
        StateSerialized() : State(TYPE) { }
        void set(int iteration, arithmetic::GWNum& X) override { TaskState::set(iteration); _serialized_value = X; }
        void to_GWNum(arithmetic::GWNum& X) override { X = _serialized_value; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_serialized_value); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_serialized_value); }

    private:
        arithmetic::SerializedGWNum _serialized_value;
    };
    class StateValue : public State
    {
    public:
        static const char TYPE = 1;
        StateValue() : State(TYPE) { }
        void set(int iteration, arithmetic::GWNum& X) override { TaskState::set(iteration); _giant_value = X; if (X.arithmetic().state().need_mod()) X.arithmetic().state().mod(_giant_value, _giant_value); }
        void to_GWNum(arithmetic::GWNum& X) override { X = _giant_value; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_giant_value); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_giant_value); }

        template<class T>
        StateValue(int iteration, T&& X) : State(TYPE) { TaskState::set(iteration); _giant_value = std::forward<T>(X); }
        template<class T>
        void set(int iteration, T&& X) { TaskState::set(iteration); _giant_value = std::forward<T>(X); }
        arithmetic::Giant& value() { return _giant_value; }

    private:
        arithmetic::Giant _giant_value;
    };

public:
    BaseExp()
    {
    }
    virtual ~BaseExp() { }

    virtual State* state() { return static_cast<State*>(Task::state()); }
    virtual arithmetic::Giant* result() { StateValue* value; if (state() == nullptr || state()->iteration() != iterations() || (value = dynamic_cast<StateValue*>(state())) == nullptr) return nullptr; return &value->value(); }
    template<class... Args>
    void commit_execute(int iteration, Args&&... args)
    {
        if (iteration == iterations())
            Task::commit_execute<StateValue>(iteration, std::forward<Args>(args)...);
        else
            Task::commit_execute<StateSerialized>(iteration, std::forward<Args>(args)...);
    }

    virtual double cost() { return _exp.bitlen(); }
    
    bool smooth() { return _smooth; }
    arithmetic::Giant& b() { return _smooth ? _exp : *(arithmetic::Giant*)nullptr; }
    arithmetic::Giant& exp() { return _exp; }
    arithmetic::Giant& tail() { return _tail; }
    arithmetic::Giant& X0() { return !_smooth ? _X0 : *(arithmetic::Giant*)nullptr; }
    uint32_t x0() { return !_smooth ? _x0 : 0; }

protected:
    bool _smooth;
    arithmetic::Giant _exp;
    arithmetic::Giant _tail;
    arithmetic::Giant _X0;
    uint32_t _x0 = 0;
};

class CarefulExp : public BaseExp
{
public:
    template<class T>
    CarefulExp(T&& exp) : BaseExp()
    {
        _smooth = false;
        _exp = std::forward<T>(exp);
    }

    void init_small(InputNum* input, arithmetic::GWState* gwstate, Logging* logging, uint32_t x0)
    {
        _x0 = x0;
        init(input, gwstate, logging);
    }
    template<class T>
    void init_small(InputNum* input, arithmetic::GWState* gwstate, Logging* logging, uint32_t x0, T&& tail)
    {
        _x0 = x0;
        _tail = std::forward<T>(tail);
        init(input, gwstate, logging);
    }
    template<class T>
    void init_giant(InputNum* input, arithmetic::GWState* gwstate, Logging* logging, T&& X0)
    {
        _X0 = std::forward<T>(X0);
        init(input, gwstate, logging);
    }
    template<class T, class S>
    void init_giant(InputNum* input, arithmetic::GWState* gwstate, Logging* logging, T&& X0, S&& tail)
    {
        _X0 = std::forward<T>(X0);
        _tail = std::forward<T>(tail);
        init(input, gwstate, logging);
    }

    double cost() override { return _exp.bitlen()*1.5; }

protected:
    void init(InputNum* input, arithmetic::GWState* gwstate, Logging* logging);
    void execute() override;
    void setup() override { }
    void release() override { }
};

class MultipointExp : public BaseExp
{
public:
    struct Point
    {
        Point(int pos_, bool check_ = true, bool value_ = true) : pos(pos_), check(check_), value(value_) {}

        int pos;
        bool check;
        bool value;
    };

public:
    template<class T>
    MultipointExp(T&& exp, bool smooth, const std::vector<Point>& points, std::function<bool(int, State*, Logging&)> on_point) : BaseExp(), _points(points), _on_point(on_point)
    {
        _exp = std::forward<T>(exp);
        _smooth = smooth;
    }

    void init_smooth(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging)
    {
        GWASSERT(smooth());
        init(input, gwstate, file, logging);
    }
    template<class T>
    void init_smooth(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& tail)
    {
        GWASSERT(smooth());
        _tail = std::forward<T>(tail);
        init(input, gwstate, file, logging);
    }
    virtual void init_state(State* state);

    void init_small(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, uint32_t x0)
    {
        GWASSERT(!smooth());
        _x0 = x0;
        init(input, gwstate, file, logging);
    }
    template<class T>
    void init_small(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, uint32_t x0, T&& tail)
    {
        GWASSERT(!smooth());
        _x0 = x0;
        _tail = std::forward<T>(tail);
        init(input, gwstate, file, logging);
    }
    template<class T>
    void init_giant(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& X0)
    {
        GWASSERT(!smooth());
        _X0 = std::forward<T>(X0);
        init(input, gwstate, file, logging);
    }
    template<class T, class S>
    void init_giant(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& X0, S&& tail)
    {
        GWASSERT(!smooth());
        _X0 = std::forward<T>(X0);
        _tail = std::forward<T>(tail);
        init(input, gwstate, file, logging);
    }

    double cost() override;
    int _W = 5;
    int _max_size = -1;
    std::vector<Point>& points() { return _points; }

protected:
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging);
    void setup() override;
    void release() override;
    void execute() override;

    int slide_init(int len);
    void slide(const arithmetic::Giant& exp, int start, int end, bool commit, int W = 0);
    void sliding_window(const arithmetic::Giant& exp);

    arithmetic::GWNum& X() { return *_X; }

protected:
    std::vector<Point> _points;
    std::function<bool(int, State*, Logging&)> _on_point;

    std::unique_ptr<arithmetic::GWNum> _X;
    std::vector<arithmetic::GWNum> _U;
};

class SmoothExp : public MultipointExp
{
public:
    SmoothExp(arithmetic::Giant& b, int n) : MultipointExp(b, true, std::vector<Point>(), nullptr)
    {
        _points.emplace_back(n);
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging)
    {
        init_smooth(input, gwstate, file, logging);
    }
    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& tail)
    {
        init_smooth(input, gwstate, file, logging, std::forward<T>(tail));
    }

private:
    using MultipointExp::init_giant;
    using MultipointExp::init_small;
    using MultipointExp::init_smooth;
};

class FastExp : public MultipointExp
{
public:
    template<class T>
    FastExp(T&& exp) : MultipointExp(std::forward<T>(exp), false, std::vector<Point>(), nullptr)
    {
        _points.emplace_back(_exp.bitlen() - 1);
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, uint32_t x0)
    {
        init_small(input, gwstate, file, logging, x0);
    }
    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, uint32_t x0, T&& tail)
    {
        init_small(input, gwstate, file, logging, x0, std::forward<T>(tail));
    }

private:
    using MultipointExp::init_giant;
    using MultipointExp::init_small;
    using MultipointExp::init_smooth;
};

class SlidingWindowExp : public MultipointExp
{
public:
    template<class T>
    SlidingWindowExp(T&& exp) : MultipointExp(std::forward<T>(exp), false, std::vector<Point>(), nullptr)
    {
        _points.emplace_back(_exp.bitlen() - 1);
    }

    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& X0)
    {
        init_giant(input, gwstate, file, logging, std::forward<T>(X0));
    }
    template<class T, class S>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& X0, S&& tail)
    {
        init_giant(input, gwstate, file, logging, std::forward<T>(X0), std::forward<T>(tail));
    }

private:
    using MultipointExp::init_giant;
    using MultipointExp::init_small;
    using MultipointExp::init_smooth;
};
 
class StrongCheckMultipointExp : public MultipointExp
{
public:
    class StrongCheckState : public TaskState
    {
    public:
        static const int TYPE = 2;
        StrongCheckState() : TaskState(TYPE) { }
        void set(int iteration, int recovery, arithmetic::GWNum& X, arithmetic::GWNum& D) { TaskState::set(iteration); _recovery = recovery; _X = X; _D = D; }
        int recovery() { return _recovery; }
        arithmetic::SerializedGWNum& X() { return _X; }
        arithmetic::SerializedGWNum& D() { return _D; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_recovery) && reader.read(_X) && reader.read(_D); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_recovery); writer.write(_X); writer.write(_D); }

    private:
        int _recovery;
        arithmetic::SerializedGWNum _X;
        arithmetic::SerializedGWNum _D;
    };

public:
    template<class T>
    StrongCheckMultipointExp(T&& exp, bool smooth, const std::vector<Point>& points, int L, int L2, std::function<bool(int, State*, Logging&)> on_point) : MultipointExp(std::forward<T>(exp), smooth, points, on_point), _L(L), _L2(L2)
    {
    }

    void init_smooth(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging)
    {
        GWASSERT(smooth());
        _tail.arithmetic().free(_tail);
        init(input, gwstate, file, file_recovery, logging);
    }
    template<class T>
    void init_smooth(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& tail)
    {
        GWASSERT(smooth());
        _tail = std::forward<T>(tail);
        _tail_inv.arithmetic().free(_tail_inv);
        init(input, gwstate, file, file_recovery, logging);
    }
    template<class T, class S>
    void init_smooth(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& tail, S&& tail_inv)
    {
        GWASSERT(smooth());
        _tail = std::forward<T>(tail);
        _tail_inv = std::forward<S>(tail_inv);
        init(input, gwstate, file, file_recovery, logging);
    }
    void init_small(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, uint32_t x0)
    {
        GWASSERT(!smooth());
        _x0 = x0;
        _tail.arithmetic().free(_tail);
        init(input, gwstate, file, file_recovery, logging);
    }
    template<class T>
    void init_small(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, uint32_t x0, T&& tail)
    {
        GWASSERT(!smooth());
        _x0 = x0;
        _tail = std::forward<T>(tail);
        _tail_inv.arithmetic().free(_tail_inv);
        init(input, gwstate, file, file_recovery, logging);
    }
    template<class T, class S>
    void init_small(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, uint32_t x0, T&& tail, S&& tail_inv)
    {
        GWASSERT(!smooth());
        _x0 = x0;
        _tail = std::forward<T>(tail);
        _tail_inv = std::forward<S>(tail_inv);
        init(input, gwstate, file, file_recovery, logging);
    }
    template<class T>
    void init_giant(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& X0)
    {
        GWASSERT(!smooth());
        _X0 = std::forward<T>(X0);
        _tail.arithmetic().free(_tail);
        init(input, gwstate, file, file_recovery, logging);
    }
    template<class T, class S>
    void init_giant(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& X0, S&& tail)
    {
        GWASSERT(!smooth());
        _X0 = std::forward<T>(X0);
        _tail = std::forward<S>(tail);
        _tail_inv.arithmetic().free(_tail_inv);
        init(input, gwstate, file, file_recovery, logging);
    }
    template<class T, class S, class R>
    void init_giant(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& X0, S&& tail, R&& tail_inv)
    {
        GWASSERT(!smooth());
        _X0 = std::forward<T>(X0);
        _tail = std::forward<S>(tail);
        _tail_inv = std::forward<R>(tail_inv);
        init(input, gwstate, file, file_recovery, logging);
    }

    virtual void init_state(State* state) override;

    State* state() override { return static_cast<State*>(_state_recovery.get()); }
    StrongCheckState* state_check() { return dynamic_cast<StrongCheckState*>(Task::state()); }
    arithmetic::GWNum& R() { return *_R; }
    arithmetic::GWNum& D() { return *_D; }
    int _L;
    int _L2;
    double cost() override;
    static void Gerbicz_params(int iters, double log2b, int& L, int &L2);

protected:
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging);
    void write_state() override;
    void setup() override;
    void release() override;
    void execute() override;

protected:
    File* _file_recovery = nullptr;
    bool _file_recovery_empty = true;
    std::unique_ptr<TaskState> _state_recovery;
    std::unique_ptr<TaskState> _tmp_state_recovery;
    int _recovery_op = 0;
    arithmetic::Giant _tail_inv;

    std::unique_ptr<arithmetic::GWNum> _R;
    std::unique_ptr<arithmetic::GWNum> _D;
};

class GerbiczCheckExp : public StrongCheckMultipointExp
{
public:
    GerbiczCheckExp(arithmetic::Giant& b, int n, int count, std::function<bool(int, State*, Logging&)> on_point = nullptr, int L = 0) : StrongCheckMultipointExp(b, true, std::vector<Point>(), 0, 0, on_point)
    {
        if (n < count)
        {
            _L = 1;
            _L2 = 1;
        }
        else if (L == 0)
            StrongCheckMultipointExp::Gerbicz_params(n/count, log2(b), _L, _L2);
        else
        {
            _L = L;
            _L2 = n/count;
            _L2 -= _L2%L;
        }
        for (int i = 0; i <= count && _L2*i <= n; i++)
            _points.emplace_back(_L2*i, true, _L2*i == n);
        if (_points.back().pos != n)
            _points.emplace_back(n);
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging)
    {
        init_smooth(input, gwstate, file, file_recovery, logging);
    }
    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& tail)
    {
        init_smooth(input, gwstate, file, file_recovery, logging, std::forward<T>(tail));
    }
    template<class T, class S>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& tail, S&& tail_inv)
    {
        init_smooth(input, gwstate, file, file_recovery, logging, std::forward<T>(tail), std::forward<S>(tail_inv));
    }

private:
    using StrongCheckMultipointExp::init_giant;
    using StrongCheckMultipointExp::init_small;
    using StrongCheckMultipointExp::init_smooth;
};

class LiCheckExp : public StrongCheckMultipointExp
{
public:
    template<class T>
    LiCheckExp(T&& exp, int count, int L = 0) : StrongCheckMultipointExp(std::forward<T>(exp), false, std::vector<Point>(), 0, 0, nullptr)
    {
        int n = _exp.bitlen() - 1;
        if (n < count)
        {
            _L = 1;
            _L2 = 1;
        }
        else if (L == 0)
            StrongCheckMultipointExp::Gerbicz_params(n/count, 1.0, _L, _L2);
        else
        {
            _L = L;
            _L2 = n/count;
            _L2 -= _L2%L;
        }
        for (int i = 0; i <= count && _L2*i <= n; i++)
            _points.emplace_back(_L2*i, true, _L2*i == n);
        if (_points.back().pos != n)
            _points.emplace_back(n);
    }

    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& X0)
    {
        init_giant(input, gwstate, file, file_recovery, logging, std::forward<T>(X0));
    }
    template<class T, class S, class R>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& X0, S&& tail, R&& tail_inv)
    {
        init_giant(input, gwstate, file, file_recovery, logging, std::forward<T>(X0), std::forward<S>(tail), std::forward<R>(tail_inv));
    }

protected:
    using StrongCheckMultipointExp::init_small;
private:
    using StrongCheckMultipointExp::init_giant;
    using StrongCheckMultipointExp::init_smooth;
};

class FastLiCheckExp : public LiCheckExp
{
public:
    template<class T>
    FastLiCheckExp(T&& exp, int count, int L = 0) : LiCheckExp(std::forward<T>(exp), count, L)
    {
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, uint32_t x0)
    {
        init_small(input, gwstate, file, file_recovery, logging, x0);
    }
    template<class T, class S>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, uint32_t x0, T&& tail, S&& tail_inv)
    {
        init_small(input, gwstate, file, file_recovery, logging, x0, std::forward<T>(tail), std::forward<S>(tail_inv));
    }

private:
    using LiCheckExp::init_small;
};

template<class IT>
class Product : public InputTask
{
public:
    Product(IT first, IT last) : InputTask(), _first(first), _last(last)
    {
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, Logging* logging)
    {
        InputTask::init(input, gwstate, nullptr, nullptr, logging, (int)(_last - _first));
    }

    arithmetic::Giant& result() { return static_cast<BaseExp::StateValue*>(state())->value(); }

protected:
    void setup() override { }
    void release() override { }
    void execute() override
    {
        int i;
        arithmetic::GWNum P(gw());
        arithmetic::GWNum X(gw());
        if (state() == nullptr)
        {
            P = *_first;
            i = 1;
        }
        else
        {
            i = state()->iteration();
            P = static_cast<BaseExp::StateValue*>(state())->value();
        }
        for (IT it = _first + i; it != _last; it++, i++, commit_execute<BaseExp::StateValue>(i, P))
        {
            X = *it;
            gw().carefully().mul(X, P, P, 0);
        }

        if (_gwstate->need_mod())
            _gwstate->mod(result(), result());
        done();
    }

protected:
    IT _first;
    IT _last;
};
