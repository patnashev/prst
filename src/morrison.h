#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "params.h"
#include "lucas.h"

class LucasMul : public InputTask
{
public:
    class State : public TaskState
    {
    public:
        static const char TYPE = 7;
        State() : TaskState(TYPE) { }
        template<class T>
        State(int iteration, int index, T&& V, bool parity) : TaskState(TYPE) { TaskState::set(iteration); _index = index; _V = std::forward<T>(V); _parity = parity; }
        template<class T>
        void set(int iteration, int index, T&& V, bool parity) { TaskState::set(iteration); _index = index; _V = std::forward<T>(V); _parity = parity; }
        int index() { return _index; }
        arithmetic::Giant& V() { return _V; }
        bool parity() { return _parity; }
        bool read(Reader& reader) override { int parity = 0; bool res = TaskState::read(reader) && reader.read(_index) && reader.read(_V) && reader.read(parity); _parity = parity == 1; return res; }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_index); writer.write(_V); writer.write(_parity ? 1 : 0); }

    private:
        int _index = 0;
        arithmetic::Giant _V;
        bool _parity;
    };

public:
    LucasMul(bool negativeQ) : InputTask(), _negativeQ(negativeQ) { }

    template<class T>
    void mul_giant(T&& b, int n)
    {
        _giants.emplace_back(std::forward<T>(b), n);
        _progress.reset();
    }
    int mul_prime(int prime, int n, int index = 0);

    std::vector<std::pair<arithmetic::Giant, int>>& giants() { return _giants; }
    std::vector<std::tuple<int, int, int>>& primes() { return _primes; }

    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* checkpoint, Logging* logging, T&& P)
    {
        init(input, gwstate, checkpoint, logging);
        if (state() == nullptr)
            init_state(new State(0, 0, std::forward<T>(P), true));
    }
    void init(InputNum* input, arithmetic::GWState* gwstate, File* checkpoint, Logging* logging);
    void init_state(State* state);

    virtual State* state() { return static_cast<State*>(Task::state()); }
    double cost();
    double progress() override;

protected:
    void setup() override;
    void release() override;
    void execute() override;

    void progress_init();

protected:
    bool _negativeQ;
    std::vector<std::pair<arithmetic::Giant, int>> _giants;
    std::vector<std::tuple<int, int, int>> _primes;
    std::unique_ptr<Progress> _progress;
};

class Morrison
{
public:
    Morrison(InputNum& input, Params& params, Logging& logging);

    void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging);

    bool success() { return _success; }
    bool prime() { return _prime; }
    std::string& res64() { return _res64; }

protected:
    class FactorTask
    {
    public:
        FactorTask(int i) : index(i) { }
        int index;
        std::unique_ptr<LucasMul> taskFactor;
        std::unique_ptr<LucasMul> taskCheck;
    };

protected:
    bool _all_factors = false;
    std::unique_ptr<LucasMul> _task;
    std::unique_ptr<LucasMul> _taskCheck;
    std::vector<FactorTask> _factor_tasks;
    int _P;
    bool _negQ;

    bool _success = false;
    bool _prime = false;
    std::string _res64;
};

