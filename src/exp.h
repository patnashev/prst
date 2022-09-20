#pragma once

#include "arithmetic.h"
#include "group.h"
#include "integer.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"

class BaseExp : public Task
{
public:
    class State : public TaskState
    {
    public:
        State() : TaskState(1) { }
        State(int iteration, const arithmetic::Giant& X) : TaskState(1) { TaskState::set(iteration); _X = X; }
        State(int iteration, arithmetic::Giant&& X) : TaskState(1) { TaskState::set(iteration); _X = std::move(X); }
        void set(int iteration, arithmetic::GWNum& X) { TaskState::set(iteration); _X = X; }
        arithmetic::Giant& X() { return _X; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_X); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_X); }

    private:
        arithmetic::Giant _X;
    };

public:
    BaseExp() : Task(true)
    {
    }
    virtual ~BaseExp() { }

    virtual State* state() { return static_cast<State*>(Task::state()); }
    double timer() { return _timer; }
    int transforms() { return _transforms; }
    void set_error_check(bool near, bool check) { _error_check_near = near; _error_check_forced = check; }

protected:
    void setup() override { }
    void release() override { }
    void reinit_gwstate() override;
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, TaskState* state, Logging* logging, int iterations);
    void done();

protected:
    InputNum* _input = nullptr;
    double _timer = 0;
    int _transforms = 0;
    bool _error_check_near = true;
    bool _error_check_forced = false;
};

class FastExp : public BaseExp
{
public:
    FastExp(const arithmetic::Giant& exp) : BaseExp(), _exp(exp) { }
    FastExp(arithmetic::Giant&& exp) : BaseExp(), _exp(std::move(exp)) { }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, uint32_t x0);
    arithmetic::Giant& exp() { return _exp; }

protected:
    void execute() override;

protected:
    arithmetic::Giant _exp;
    uint32_t _x0;
};

class SlowExp : public BaseExp
{
public:
    SlowExp(const arithmetic::Giant& exp) : BaseExp(), _exp(exp) { }
    SlowExp(arithmetic::Giant&& exp) : BaseExp(), _exp(std::move(exp)) { }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, const arithmetic::Giant& X0);
    arithmetic::Giant& exp() { return _exp; }

protected:
    void execute() override;

protected:
    arithmetic::Giant _exp;
    arithmetic::Giant _X0;
};

class MultipointExp : public BaseExp
{
public:
    MultipointExp(arithmetic::Giant& b, const std::vector<int>& points, std::function<void(int)> on_point) : BaseExp(), _b(b), _points(points), _on_point(on_point)
    {
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging);
    virtual void init_state(State* state);

    arithmetic::GWNum& X() { return *_X; }
    int _W = 5;
    int _max_size = -1;

protected:
    void release() override;
    void execute() override;
    void sliding_window(const arithmetic::Giant& exp);

protected:
    arithmetic::Giant& _b;
    const std::vector<int>& _points;
    std::function<void(int)> _on_point;

    arithmetic::Giant _X0;
    std::unique_ptr<arithmetic::GWNum> _X;
    std::vector<arithmetic::GWNum> _U;
};

class GerbiczCheckMultipointExp : public MultipointExp
{
public:
    static int CHECKS_PER_POINT;

public:
    class GerbiczCheckState : public TaskState
    {
    public:
        GerbiczCheckState() : TaskState(2) { }
        void set(int iteration, arithmetic::GWNum& X, arithmetic::GWNum& D) { TaskState::set(iteration); _X = X; _D = D; }
        arithmetic::Giant& X() { return _X; }
        arithmetic::Giant& D() { return _D; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_X) && reader.read(_D); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_X); writer.write(_D); }

    private:
        arithmetic::Giant _X;
        arithmetic::Giant _D;
    };

public:
    GerbiczCheckMultipointExp(arithmetic::Giant& b, const std::vector<int>& points, std::function<void(int)> on_point) : MultipointExp(b, points, on_point)
    {
        int iters = points[0];
        for (int i = 1; i < points.size() - 1; i++)
            if (iters > points[i] - points[i - 1])
                iters = points[i] - points[i - 1];
        Gerbicz_params(iters/CHECKS_PER_POINT, log2(b), _L, _L2);
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging);
    void init_state(State* state) override;

    State* state() override { return _state_recovery.get(); }
    GerbiczCheckState* state_check() { return dynamic_cast<GerbiczCheckState*>(Task::state()); }
    arithmetic::GWNum& R() { return *_R; }
    arithmetic::GWNum& D() { return *_D; }
    int _L;
    int _L2;
    double cost();
    static void Gerbicz_params(int iters, double log2b, int& L, int &L2);

protected:
    void write_state() override;
    void release() override;
    void setup() override;
    void execute() override;

protected:
    File* _file_recovery = nullptr;
    std::unique_ptr<State> _state_recovery;
    int _recovery_op = 0;

    std::unique_ptr<arithmetic::GWNum> _R;
    std::unique_ptr<arithmetic::GWNum> _D;
};

class GerbiczCheckPoints
{
public:
    GerbiczCheckPoints(double log2b, int n, int count)
    {
        int L, L2;
        GerbiczCheckMultipointExp::Gerbicz_params(n/count, log2b, L, L2);
        for (int i = 1; i <= count; i++)
            _recovery_points.push_back(L2*i);
        if (_recovery_points.back() != n)
            _recovery_points.push_back(n);
    }
    GerbiczCheckPoints(int n, int count, int L)
    {
        int L2 = n/count;
        L2 -= L2%L;
        for (int i = 1; i <= count; i++)
            _recovery_points.push_back(L2*i);
        if (_recovery_points.back() != n)
            _recovery_points.push_back(n);
    }

protected:
    std::vector<int> _recovery_points;
};

class GerbiczCheckExp : private GerbiczCheckPoints, public GerbiczCheckMultipointExp
{
public:
    GerbiczCheckExp(arithmetic::Giant& b, int n, int count) : GerbiczCheckPoints(log2(b), n, count), GerbiczCheckMultipointExp(b, _recovery_points, nullptr)
    {
    }
    GerbiczCheckExp(arithmetic::Giant& b, int n, int count, int L) : GerbiczCheckPoints(n, count, L), GerbiczCheckMultipointExp(b, _recovery_points, nullptr)
    {
        _L = L;
        _L2 = n/count;
        _L2 -= _L2%L;
    }
};
