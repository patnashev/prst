#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "options.h"
#include "exp.h"

class Proof
{
public:
    static const int NO_OP = 0;
    static const int SAVE = 1;
    static const int BUILD = 2;
    static const int CERT = 3;
    static const int ROOT = 4;

public:
    class Product : public TaskState
    {
    public:
        static const char TYPE = 3;
        Product() : TaskState(TYPE) { }
        void set(int power, arithmetic::GWNum& X) { TaskState::set(power); _X = X; if (X.arithmetic().state().need_mod()) X.arithmetic().state().mod(_X, _X); }
        int depth() { return _iteration; }
        arithmetic::Giant& value() { return _X; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_X); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_X); }

    private:
        arithmetic::Giant _X;
    };
    class Certificate : public TaskState
    {
    public:
        static const char TYPE = 4;
        Certificate() : TaskState(TYPE) { }
        void set(int power, arithmetic::GWNum& X) { TaskState::set(power); _X = X; _a_power = 0; _a_base = 0; if (X.arithmetic().state().need_mod()) X.arithmetic().state().mod(_X, _X); }
        void set(int power, arithmetic::GWNum& X, arithmetic::Giant& a_power, arithmetic::GWNum& A) { TaskState::set(power); _X = X; _a_power = a_power; _a_base = A; if (X.arithmetic().state().need_mod()) { X.arithmetic().state().mod(_X, _X); A.arithmetic().state().mod(_a_base, _a_base); } }
        int power() { return _iteration; }
        arithmetic::Giant& X() { return _X; }
        arithmetic::Giant& a_power() { return _a_power; }
        arithmetic::Giant& a_base() { return _a_base; }
        bool read(Reader& reader) override { _a_power = 0; _a_base = 0; return TaskState::read(reader) && reader.read(_X) && ((reader.read(_a_power) && _a_power != 0 && reader.read(_a_base)) || true); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_X); writer.write(_a_power); writer.write(_a_base); }

    private:
        arithmetic::Giant _X;
        arithmetic::Giant _a_power;
        arithmetic::Giant _a_base;
    };
    class State : public TaskState
    {
    public:
        static const char TYPE = 6;
        State() : TaskState(TYPE) { }
        template<class T>
        State(int iteration, T&& Y) : TaskState(TYPE) { TaskState::set(iteration); _Y = std::forward<T>(Y); }
        template<class T, class R, class S>
        State(int iteration, S&& X, T&& Y, R&& exp) : TaskState(TYPE) { TaskState::set(iteration); _X = std::forward<S>(X); _Y = std::forward<T>(Y); _exp = std::forward<R>(exp); }
        void set(int iteration, arithmetic::GWNum& Y, const std::vector<arithmetic::Giant>& h) { TaskState::set(iteration); _Y = Y; _h = h; if (Y.arithmetic().state().need_mod()) Y.arithmetic().state().mod(_Y, _Y); }
        void set(int iteration, arithmetic::GWNum& X, arithmetic::GWNum& Y, arithmetic::Giant& exp, const std::vector<arithmetic::Giant>& h) { TaskState::set(iteration); _X = X; _Y = Y; _exp = exp; _h = h; if (Y.arithmetic().state().need_mod()) Y.arithmetic().state().mod(_Y, _Y); }
        arithmetic::SerializedGWNum& X() { return _X; }
        arithmetic::Giant& Y() { return _Y; }
        arithmetic::Giant& exp() { return _exp; }
        std::vector<arithmetic::Giant>& h() { return _h; }

    private:
        arithmetic::SerializedGWNum _X;
        arithmetic::Giant _Y;
        arithmetic::Giant _exp;
        std::vector<arithmetic::Giant> _h;
    };


public:
    Proof(int op, int count, InputNum& input, Options& options, File& file_cert, Logging& logging);

    void calc_points(int iterations, bool smooth, InputNum& input, Options& options, Logging& logging);
    void init_files(File* file_point, File* file_product, File* file_cert);
    void init_state(MultipointExp* task, arithmetic::GWState& gwstate, InputNum& input, Logging& logging, int a);
    void read_point(int index, TaskState& state, Logging& logging);
    BaseExp::State* read_point(int index, BaseExp::StateValue* state_value, BaseExp::StateSerialized* state_serialized, Logging& logging);
    void read_product(int index, TaskState& state, Logging& logging);
    bool on_point(int index, BaseExp::State* state, Logging& logging);
    void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging);
    void run(InputNum& input, arithmetic::GWState& gwstate, Logging& logging, arithmetic::Giant* X);
    double cost();

    int op() { return _op; }
    int count() { return _count; }
    bool Li() { return _Li; }
    int depth() { int t = 0; while ((1 << t) < _count) t++; return t; }
    std::vector<MultipointExp::Point>& points() { return _points; }
    int M() { return _M; }
    void set_cache_points(bool value) { _cache_points = value; }
    InputTask* task() { return _task.get(); }
    CarefulExp* taskRoot() { return _taskRoot.get(); }
    std::string& res64() { return _res64; }

    std::vector<File*>& file_points() { return _file_points; }
    std::vector<File*>& file_products() { return _file_products; }
    File* file_cert() { return _file_cert; }
    arithmetic::Giant& r_0() { return _r_0; }
    arithmetic::Giant& r_count() { return _r_count; }
    arithmetic::Giant& r_exp() { return *_r_exp; }

protected:

protected:
    int _op = 0;
    int _count = 0;
    bool _Li = false;
    std::vector<MultipointExp::Point> _points;
    int _M = 0;
    bool _cache_points = false;
    std::vector<File*> _file_points;
    std::vector<File*> _file_products;
    File* _file_cert = nullptr;
    arithmetic::Giant _r_0;
    arithmetic::Giant _r_count;
    arithmetic::Giant* _r_exp;
    std::unique_ptr<InputTask> _task;
    std::unique_ptr<MultipointExp> _taskA;
    std::unique_ptr<CarefulExp> _taskRoot;
    std::string _res64;
};

class ProofSave : public InputTask
{
public:
    ProofSave() { }

    void init(InputNum* input, arithmetic::GWState* gwstate, Logging* logging, Proof* proof);

    Proof::State* state() { return static_cast<Proof::State*>(Task::state()); }

protected:
    void setup() override { }
    void release() override { }
    void execute() override;
    void done() override;

protected:
    Proof* _proof;
};

class ProofBuild : public InputTask
{
public:
    ProofBuild(const std::string& security_seed) : _security_seed(security_seed) { }

    void init(InputNum* input, arithmetic::GWState* gwstate, Logging* logging, Proof* proof);

    Proof::State* state() { return static_cast<Proof::State*>(Task::state()); }
    bool security() { return !_security_seed.empty(); }
    std::string& raw_res64() { return _raw_res64; }

protected:
    void setup() override { }
    void release() override { }
    void execute() override;
    void done() override;

protected:
    Proof* _proof;
    std::string _security_seed;
    arithmetic::Giant _rnd_seed;
    std::string _raw_res64;
};
