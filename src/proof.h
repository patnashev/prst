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
    class Certificate : public TaskState
    {
    public:
        static const int TYPE = 3;
        Certificate() : TaskState(TYPE) { }
        void set(int power, arithmetic::GWNum& X) { TaskState::set(power); _X = X; }
        int power() { return _iteration; }
        arithmetic::Giant& X() { return _X; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_X); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_X); }

    private:
        arithmetic::Giant _X;
    };

public:
    Proof(int op, int count, InputNum& input, Params& params, Logging& logging);

    static int read_cert_power(File& file_cert);
    void calc_points(int iterations, InputNum& input, Params& params, Logging& logging);
    void init_files(File* file_point, File* file_product, File* file_cert);
    void init_state(MultipointExp* task, arithmetic::GWState& gwstate, InputNum& input, Logging& logging, int a);
    void read_point(int index, TaskState& state, Logging& logging);
    void read_product(int index, TaskState& state, Logging& logging);
    void on_point(int index, arithmetic::Giant& X);
    void run(InputNum& input, arithmetic::GWState& gwstate, File& file_cert, File& file_checkpoint, File& file_recoverypoint, Logging& logging);
    void run(InputNum& input, arithmetic::GWState& gwstate, Logging& logging, arithmetic::Giant* X);
    double cost();

    int op() { return _op; }
    int count() { return _count; }
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

protected:

protected:
    int _op = 0;
    int _count = 0;
    std::vector<int> _points;
    int _M = 0;
    bool _cache_points = false;
    std::vector<File*> _file_points;
    std::vector<File*> _file_products;
    File* _file_cert = nullptr;
    arithmetic::Giant _r_0;
    arithmetic::Giant _r_count;
    std::unique_ptr<BaseExp> _task;
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
