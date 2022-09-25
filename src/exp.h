#pragma once

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
        static const int TYPE = 1;
        State() : TaskState(TYPE) { }
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
    BaseExp()
    {
    }
    virtual ~BaseExp() { }

    virtual State* state() { return static_cast<State*>(Task::state()); }
    double timer() { return _timer; }
    int transforms() { return _transforms; }
    void set_error_check(bool near, bool check);

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
    template<class T>
    FastExp(T&& exp) : BaseExp()
    {
        _exp = std::forward<T>(exp);
    }

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
    virtual double cost();

protected:
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging);
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
        void set(int iteration, int recovery, arithmetic::GWNum& X, arithmetic::GWNum& D) { TaskState::set(iteration); _recovery = recovery; _X = X; _D = D; }
        int recovery() { return _recovery; }
        arithmetic::Giant& X() { return _X; }
        arithmetic::Giant& D() { return _D; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_recovery) && reader.read(_X) && reader.read(_D); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_recovery); writer.write(_X); writer.write(_D); }

    private:
        int _recovery;
        arithmetic::Giant _X;
        arithmetic::Giant _D;
    };

public:
    GerbiczCheckMultipointExp(arithmetic::Giant& b, const std::vector<int>& points, int points_per_check, int checks_per_point, std::function<void(int, arithmetic::Giant&)> on_point) : MultipointExp(b, points, on_point), _points_per_check(points_per_check), _checks_per_point(checks_per_point)
    {
        Gerbicz_params(points[points_per_check < points.size() ? points_per_check : points.size() - 1]/checks_per_point, log2(b), _L, _L2);
    }

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
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging);
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
        for (int i = 0; i <= count; i++)
            _recovery_points.push_back(L2*i);
        if (_recovery_points.back() != n)
            _recovery_points.push_back(n);
    }
    GerbiczCheckPoints(int n, int count, int L)
    {
        int L2 = n/count;
        L2 -= L2%L;
        for (int i = 0; i <= count; i++)
            _recovery_points.push_back(L2*i);
        if (_recovery_points.back() != n)
            _recovery_points.push_back(n);
    }

    std::vector<int>& recovery_points() { return _recovery_points; }

protected:
    std::vector<int> _recovery_points;
};

class GerbiczCheckExp : private GerbiczCheckPoints, public GerbiczCheckMultipointExp
{
public:
    GerbiczCheckExp(arithmetic::Giant& b, int n, int count) : GerbiczCheckPoints(log2(b), n, count), GerbiczCheckMultipointExp(b, _recovery_points, 1, 1, nullptr)
    {
    }
    GerbiczCheckExp(arithmetic::Giant& b, int n, int count, int L) : GerbiczCheckPoints(n, count, L), GerbiczCheckMultipointExp(b, _recovery_points, 1, 1, nullptr)
    {
        _L = L;
        _L2 = n/count;
        _L2 -= _L2%L;
    }
};
