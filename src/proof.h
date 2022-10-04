#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "params.h"
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
        void mimic_type(int type) { _type = type; }
        void set(int power, arithmetic::GWNum& X) { TaskState::set(power); _X = X; }
        int depth() { return _iteration; }
        arithmetic::Giant& X() { return _X; }
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
        void set(int power, arithmetic::GWNum& X) { TaskState::set(power); _X = X; _a_power = 0; _a_base = 0; }
        void set(int power, arithmetic::GWNum& X, arithmetic::Giant& a_power, arithmetic::GWNum& A) { TaskState::set(power); _X = X; _a_power = a_power; _a_base = A; }
        int power() { return _iteration; }
        arithmetic::Giant& X() { return _X; }
        arithmetic::Giant& a_power() { return _a_power; }
        arithmetic::Giant& a_base() { return _a_base; }
        bool read(Reader& reader) override { _a_power = 0; _a_base = 0; return TaskState::read(reader) && reader.read(_X) && ((reader.read(_a_power) && reader.read(_a_base)) || true); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_X); writer.write(_a_power); writer.write(_a_base); }

    private:
        arithmetic::Giant _X;
        arithmetic::Giant _a_power;
        arithmetic::Giant _a_base;
    };

public:
    Proof(int op, int count, InputNum& input, Params& params, File& file_cert, Logging& logging, std::optional<bool> Li = std::optional<bool>());

    void calc_points(int iterations, InputNum& input, Params& params, Logging& logging);
    void init_files(File* file_point, File* file_product, File* file_cert);
    void init_state(MultipointExp* task, arithmetic::GWState& gwstate, InputNum& input, Logging& logging, int a);
    void read_point(int index, TaskState& state, Logging& logging);
    void read_product(int index, TaskState& state, Logging& logging);
    bool on_point(int index, arithmetic::Giant& X);
    void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging);
    void run(InputNum& input, arithmetic::GWState& gwstate, Logging& logging, arithmetic::Giant* X);
    double cost();

    int op() { return _op; }
    int count() { return _count; }
    bool Li() { return _Li; }
    int depth() { int t; for (t = 0; (1 << t) < _count; t++); return t; }
    std::vector<int>& points() { return _points; }
    int M() { return _M; }
    BaseExp* task() { return _task.get(); }
    BaseExp* taskRoot() { return _taskRoot.get(); }
    std::string& res64() { return _res64; }

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
    std::vector<int> _points;
    int _points_per_check = 1;
    int _M = 0;
    bool _cache_points = false;
    std::vector<File*> _file_points;
    std::vector<File*> _file_products;
    File* _file_cert = nullptr;
    arithmetic::Giant _r_0;
    arithmetic::Giant _r_count;
    arithmetic::Giant* _r_exp;
    std::unique_ptr<BaseExp> _task;
    std::unique_ptr<BaseExp> _taskA;
    std::unique_ptr<BaseExp> _taskRoot;
    std::string _res64;
};

class ProofSave : public BaseExp
{
public:
    ProofSave(Proof& proof) : _proof(proof)
    {
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, Logging* logging);

protected:
    void execute() override;
    void read_point(int index, TaskState& state);

protected:
    Proof& _proof;
};

class ProofBuild : public BaseExp
{
public:
    ProofBuild(Proof& proof, const std::string& security_seed) : _proof(proof), _security_seed(security_seed)
    {
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, Logging* logging);

    bool security() { return !_security_seed.empty(); }
    std::string& raw_res64() { return _raw_res64; }

protected:
    void execute() override;

protected:
    Proof& _proof;
    std::string _security_seed;
    arithmetic::Giant _rnd_seed;
    std::string _raw_res64;
};
