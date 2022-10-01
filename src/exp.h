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
        State(int iteration, const arithmetic::Giant& X) : TaskState(1) { TaskState::set(iteration); _X = X; }
        State(int iteration, arithmetic::Giant&& X) : TaskState(1) { TaskState::set(iteration); _X = std::move(X); }
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
    double timer() { return _timer; }
    int transforms() { return _transforms; }
    void set_error_check(bool near, bool check);
    virtual double cost() { return 0; }

protected:
    void setup() override { }
    void release() override { }
    void reinit_gwstate() override;
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, TaskState* state, Logging* logging, int iterations);
    void done();

protected:
    double _timer = 0;
    int _transforms = 0;
    bool _error_check_near = true;
    bool _error_check_forced = false;
};

class FastExp : public BaseExp
{
public:
    FastExp(FastExp&& a)
    {
        _exp = std::move(a._exp);
    }
    template<class T>
    FastExp(T&& exp) : BaseExp()
    {
        _exp = std::forward<T>(exp);
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, uint32_t x0);
    arithmetic::Giant& exp() { return _exp; }
    double cost() override { return _exp.bitlen(); }

protected:
    void execute() override;

protected:
    arithmetic::Giant _exp;
    uint32_t _x0;
};

class SlowExp : public BaseExp
{
public:
    SlowExp(SlowExp&& a)
    {
        _exp = std::move(a._exp);
    }
    template<class T>
    SlowExp(T&& exp) : BaseExp()
    {
        _exp = std::forward<T>(exp);
    }

    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& X0)
    {
        _X0 = std::forward<T>(X0);
        init(input, gwstate, file, logging);
    }
    arithmetic::Giant& exp() { return _exp; }
    arithmetic::Giant& X0() { return _X0; }
    double cost() override { return _exp.bitlen()*1.5; }

protected:
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging);
    void execute() override;

protected:
    arithmetic::Giant _exp;
    arithmetic::Giant _X0;
};

class MultipointExp : public BaseExp
{
public:
    MultipointExp(arithmetic::Giant& b, const std::vector<int>& points, std::function<void(int, arithmetic::Giant&)> on_point) : BaseExp(), _b(b), _points(points), _on_point(on_point)
    {
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging);
    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, T&& tail)
    {
        _tail = std::forward<T>(tail);
        init(input, gwstate, file, logging);
    }
    virtual void init_state(State* state);

    arithmetic::GWNum& X() { return *_X; }
    int _W = 5;
    int _max_size = -1;
    double cost() override;

protected:
    void release() override;
    void execute() override;
    void sliding_window(const arithmetic::Giant& exp);

protected:
    arithmetic::Giant& _b;
    const std::vector<int>& _points;
    std::function<void(int, arithmetic::Giant&)> _on_point;

    arithmetic::Giant _tail;
    std::unique_ptr<arithmetic::GWNum> _X;
    std::vector<arithmetic::GWNum> _U;
};

class GerbiczCheckMultipointExp : public MultipointExp
{
public:
    class GerbiczCheckState : public TaskState
    {
    public:
        static const int TYPE = 2;
        GerbiczCheckState() : TaskState(TYPE) { }
        GerbiczCheckState(arithmetic::GWArithmetic& gw) : TaskState(TYPE) { _X.reset(new arithmetic::GWNum(gw)); _D.reset(new arithmetic::GWNum(gw)); }
        void set(int iteration, int recovery, arithmetic::GWNum& X, arithmetic::GWNum& D);
        int recovery() { return _recovery; }
        arithmetic::GWNum& X() { return *_X; }
        arithmetic::GWNum& D() { return *_D; }
        bool read(Reader& reader) override;
        void write(Writer& writer) override;

    private:
        int _recovery;
        std::unique_ptr<arithmetic::GWNum> _X;
        std::unique_ptr<arithmetic::GWNum> _D;
    };

public:
    GerbiczCheckMultipointExp(arithmetic::Giant& b, const std::vector<int>& points, int L, int L2, int points_per_check, int checks_per_point, std::function<void(int, arithmetic::Giant&)> on_point) : MultipointExp(b, points, on_point), _L(L), _L2(L2), _points_per_check(points_per_check), _checks_per_point(checks_per_point)
    {
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging);
    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, T&& tail)
    {
        _tail = std::forward<T>(tail);
        init(input, gwstate, file, file_recovery, logging);
    }
    void init_state(State* state) override;

    State* state() override { return _state_recovery.get(); }
    GerbiczCheckState* state_check() { return dynamic_cast<GerbiczCheckState*>(Task::state()); }
    arithmetic::GWNum& R() { return *_R; }
    arithmetic::GWNum& D() { return *_D; }
    int _L;
    int _L2;
    int _points_per_check;
    int _checks_per_point;
    double cost() override;
    static void Gerbicz_params(int iters, double log2b, int& L, int &L2);

protected:
    void write_state() override;
    void release() override;
    void reinit_gwstate() override;
    void setup() override;
    void execute() override;

protected:
    File* _file_recovery = nullptr;
    std::unique_ptr<State> _state_recovery;
    std::unique_ptr<State> _tmp_state_recovery;
    int _recovery_op = 0;

    std::unique_ptr<arithmetic::GWNum> _R;
    std::unique_ptr<arithmetic::GWNum> _D;
};

class GerbiczCheckPoints
{
public:
    GerbiczCheckPoints(double log2b, int n, int count)
    {
        GerbiczCheckMultipointExp::Gerbicz_params(n/count, log2b, _L, _L2);
        for (int i = 0; i <= count; i++)
            _recovery_points.push_back(_L2*i);
        if (_recovery_points.back() != n)
            _recovery_points.push_back(n);
    }
    GerbiczCheckPoints(int n, int count, int L) : _L(L)
    {
        _L2 = n/count;
        _L2 -= _L2%L;
        for (int i = 0; i <= count; i++)
            _recovery_points.push_back(_L2*i);
        if (_recovery_points.back() != n)
            _recovery_points.push_back(n);
    }

    std::vector<int>& recovery_points() { return _recovery_points; }

protected:
    int _L, _L2;
    std::vector<int> _recovery_points;
};

class GerbiczCheckExp : public GerbiczCheckPoints, public GerbiczCheckMultipointExp
{
public:
    GerbiczCheckExp(arithmetic::Giant& b, int n, int count, std::function<void(int, arithmetic::Giant&)> on_point = nullptr) : GerbiczCheckPoints(log2(b), n, count), GerbiczCheckMultipointExp(b, _recovery_points, GerbiczCheckPoints::_L, GerbiczCheckPoints::_L2, 1, 1, on_point)
    {
    }
    GerbiczCheckExp(arithmetic::Giant& b, int n, int count, int L, std::function<void(int, arithmetic::Giant&)> on_point = nullptr) : GerbiczCheckPoints(n, count, L), GerbiczCheckMultipointExp(b, _recovery_points, GerbiczCheckPoints::_L, GerbiczCheckPoints::_L2, 1, 1, on_point)
    {
    }
};
