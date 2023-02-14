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
    class State : public TaskState
    {
    public:
        static const char TYPE = 1;
        State() : TaskState(TYPE) { }
        State(int iteration, const arithmetic::Giant& X) : TaskState(TYPE) { TaskState::set(iteration); _X = X; }
        State(int iteration, arithmetic::Giant&& X) : TaskState(TYPE) { TaskState::set(iteration); _X = std::move(X); }
        template<class T>
        void set(int iteration, T& X) { TaskState::set(iteration); _X = X; }
        arithmetic::Giant& X() { return _X; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_X); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_X); }

    private:
        arithmetic::Giant _X;
    };

public:
    BaseExp()
    {
    }
    virtual ~BaseExp() { }

    virtual State* state() { return static_cast<State*>(Task::state()); }
    virtual double cost() { return _exp.bitlen(); }

    bool smooth() { return _smooth; }
    arithmetic::Giant& b() { return _smooth ? _exp : *(arithmetic::Giant*)nullptr; }
    arithmetic::Giant& exp() { return _exp; }
    arithmetic::Giant& tail() { return _tail; }
    arithmetic::Giant& X0() { return !_smooth ? _X0 : *(arithmetic::Giant*)nullptr; }
    uint32_t x0() { return !_smooth ? _x0 : 0; }

protected:
    void done() override;

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
    void init(InputNum* input, arithmetic::GWState* gwstate, Logging* logging, T&& X0)
    {
        _X0 = std::forward<T>(X0);
        init(input, gwstate, logging);
    }
    template<class T, class S>
    void init(InputNum* input, arithmetic::GWState* gwstate, Logging* logging, T&& X0, S&& tail)
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
    template<class T>
    MultipointExp(T&& exp, bool smooth, const std::vector<int>& points, std::function<bool(int, arithmetic::Giant&)> on_point) : BaseExp(), _points(points), _on_point(on_point)
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
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& X0)
    {
        GWASSERT(!smooth());
        _X0 = std::forward<T>(X0);
        init(input, gwstate, file, logging);
    }
    template<class T, class S>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& X0, S&& tail)
    {
        GWASSERT(!smooth());
        _X0 = std::forward<T>(X0);
        _tail = std::forward<T>(tail);
        init(input, gwstate, file, logging);
    }

    double cost() override;
    int _W = 5;
    int _max_size = -1;
    std::vector<int>& points() { return _points; }

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
    std::vector<int> _points;
    std::function<bool(int, arithmetic::Giant&)> _on_point;

    std::unique_ptr<arithmetic::GWNum> _X;
    std::vector<arithmetic::GWNum> _U;
};

class SmoothExp : public MultipointExp
{
public:
    SmoothExp(arithmetic::Giant& b, int n) : MultipointExp(b, true, std::vector<int>(), nullptr)
    {
        _points.push_back(n);
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging)
    {
        MultipointExp::init(input, gwstate, file, logging);
    }
    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& tail)
    {
        _tail = std::forward<T>(tail);
        MultipointExp::init(input, gwstate, file, logging);
    }
};

class FastExp : public MultipointExp
{
public:
    template<class T>
    FastExp(T&& exp) : MultipointExp(std::forward<T>(exp), false, std::vector<int>(), nullptr)
    {
        _points.push_back(_exp.bitlen() - 1);
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, uint32_t x0)
    {
        _x0 = x0;
        MultipointExp::init(input, gwstate, file, logging);
    }
    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, uint32_t x0, T&& tail)
    {
        _x0 = x0;
        _tail = std::forward<T>(tail);
        MultipointExp::init(input, gwstate, file, logging);
    }
};

class SlidingWindowExp : public MultipointExp
{
public:
    template<class T>
    SlidingWindowExp(T&& exp) : MultipointExp(std::forward<T>(exp), false, std::vector<int>(), nullptr)
    {
        _points.push_back(_exp.bitlen() - 1);
    }

    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& X0)
    {
        _X0 = std::forward<T>(X0);
        MultipointExp::init(input, gwstate, file, logging);
    }
    template<class T, class S>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& X0, S&& tail)
    {
        _X0 = std::forward<T>(X0);
        _tail = std::forward<T>(tail);
        MultipointExp::init(input, gwstate, file, logging);
    }
};
 
class StrongCheckMultipointExp : public MultipointExp
{
public:
    class StrongCheckState : public TaskState
    {
    public:
        static const int TYPE = 2;
        StrongCheckState() : TaskState(TYPE) { }
        void set(int iteration, int recovery, arithmetic::GWNum& X, arithmetic::GWNum& D);
        int recovery() { return _recovery; }
        arithmetic::Giant& X() { return _X; }
        arithmetic::Giant& D() { return _D; }
        std::unique_ptr<arithmetic::GWNum>& gwX() { return _gwX; }
        std::unique_ptr<arithmetic::GWNum>& gwD() { return _gwD; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_recovery) && reader.read(_X) && reader.read(_D); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_recovery); writer.write(_X); writer.write(_D); }

    private:
        int _recovery;
        arithmetic::Giant _X;
        arithmetic::Giant _D;
        std::unique_ptr<arithmetic::GWNum> _gwX;
        std::unique_ptr<arithmetic::GWNum> _gwD;
    };

public:
    template<class T>
    StrongCheckMultipointExp(T&& exp, bool smooth, const std::vector<int>& points, int L, int L2, std::function<bool(int, arithmetic::Giant&)> on_point) : MultipointExp(std::forward<T>(exp), smooth, points, on_point), _L(L), _L2(L2)
    {
    }

    void init_smooth(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging)
    {
        GWASSERT(smooth());
        init(input, gwstate, file, file_recovery, logging);
    }
    template<class T>
    void init_smooth(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& tail)
    {
        GWASSERT(smooth());
        _tail = std::forward<T>(tail);
        init(input, gwstate, file, file_recovery, logging);
    }
    virtual void init_state(State* state);

    void init_small(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, uint32_t x0)
    {
        GWASSERT(!smooth());
        _x0 = x0;
        init(input, gwstate, file, file_recovery, logging);
    }
    template<class T>
    void init_small(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, uint32_t x0, T&& tail)
    {
        GWASSERT(!smooth());
        _x0 = x0;
        _tail = std::forward<T>(tail);
        init(input, gwstate, file, file_recovery, logging);
    }
    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& X0)
    {
        GWASSERT(!smooth());
        _X0 = std::forward<T>(X0);
        init(input, gwstate, file, file_recovery, logging);
    }
    template<class T, class S>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& X0, S&& tail)
    {
        GWASSERT(!smooth());
        _X0 = std::forward<T>(X0);
        _tail = std::forward<T>(tail);
        init(input, gwstate, file, file_recovery, logging);
    }

    State* state() override { return _state_recovery.get(); }
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
    std::unique_ptr<State> _state_recovery;
    std::unique_ptr<State> _tmp_state_recovery;
    int _recovery_op = 0;

    std::unique_ptr<arithmetic::GWNum> _R;
    std::unique_ptr<arithmetic::GWNum> _D;
};

class GerbiczCheckExp : public StrongCheckMultipointExp
{
public:
    GerbiczCheckExp(arithmetic::Giant& b, int n, int count, std::function<bool(int, arithmetic::Giant&)> on_point = nullptr, int L = 0) : StrongCheckMultipointExp(b, true, std::vector<int>(), 0, 0, on_point)
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
            _points.push_back(_L2*i);
        if (_points.back() != n)
            _points.push_back(n);
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging)
    {
        StrongCheckMultipointExp::init(input, gwstate, file, file_recovery, logging);
    }
    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& tail)
    {
        _tail = std::forward<T>(tail);
        StrongCheckMultipointExp::init(input, gwstate, file, file_recovery, logging);
    }
};

class LiCheckExp : public StrongCheckMultipointExp
{
public:
    template<class T>
    LiCheckExp(T&& exp, int count, int L = 0) : StrongCheckMultipointExp(std::forward<T>(exp), false, std::vector<int>(), 0, 0, nullptr)
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
            _points.push_back(_L2*i);
        if (_points.back() != n)
            _points.push_back(n);
    }

    using StrongCheckMultipointExp::init_small;
    using StrongCheckMultipointExp::init;
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
        _x0 = x0;
        StrongCheckMultipointExp::init(input, gwstate, file, file_recovery, logging);
    }
    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, uint32_t x0, T&& tail)
    {
        _x0 = x0;
        _tail = std::forward<T>(tail);
        StrongCheckMultipointExp::init(input, gwstate, file, file_recovery, logging);
    }
};
