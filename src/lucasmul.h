#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "lucas.h"

class LucasVMul : public InputTask
{
public:

    virtual double cost() = 0;

    virtual arithmetic::Giant* result() = 0;

protected:
    void done() override;
};

class LucasVMulFast : public LucasVMul
{
public:
    class State : public TaskState
    {
    public:
        static const char TYPE = 9;
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
    LucasVMulFast(bool carefully = false) : LucasVMul(), _carefully(carefully) { }

    template<class T>
    void mul_giant(T&& b, int n)
    {
        _giants.emplace_back(std::forward<T>(b), n);
        _progress.reset();
    }
    int mul_prime(int prime, int n, int index = 0);

    std::vector<std::pair<arithmetic::Giant, int>>& giants() { return _giants; }
    std::vector<std::tuple<int, int, int>>& primes() { return _primes; }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* checkpoint, Logging* logging, bool negativeQ)
    {
        _negativeQ = negativeQ;
        init(input, gwstate, checkpoint, logging);
    }
    template<class T>
    void init(InputNum* input, arithmetic::GWState* gwstate, File* checkpoint, Logging* logging, T&& P, bool negativeQ)
    {
        _negativeQ = negativeQ;
        init(input, gwstate, checkpoint, logging);
        if (state() == nullptr)
            init_state(new State(0, 0, std::forward<T>(P), true));
    }
    void init_state(State* state);

    virtual State* state() { return static_cast<State*>(Task::state()); }
    double cost() override;
    double progress() override;
    arithmetic::Giant* result() override { if (state() == nullptr || state()->index() != _giants.size() + _primes.size()) return nullptr; return &state()->V(); }

protected:
    void init(InputNum* input, arithmetic::GWState* gwstate, File* checkpoint, Logging* logging);
    void setup() override;
    void release() override;
    void execute() override;

    void progress_init();

protected:
    bool _carefully;
    bool _negativeQ;
    std::vector<std::pair<arithmetic::Giant, int>> _giants;
    std::vector<std::tuple<int, int, int>> _primes;
    std::unique_ptr<Progress> _progress;
};

class LucasUVMulFast : public LucasVMul
{
public:
    class State : public TaskState
    {
    public:
        static const char TYPE = 10;
        State() : TaskState(TYPE) { }
        template<class T0, class T1>
        State(int iteration, T0&& Vn, T1&& Vn1, bool parity) : TaskState(TYPE) { TaskState::set(iteration); _Vn = std::forward<T0>(Vn); _Vn1 = std::forward<T1>(Vn1); _parity = parity; }
        void set(int iteration, const arithmetic::LucasV& Vn, const arithmetic::LucasV& Vn1) { TaskState::set(iteration); _Vn = Vn.V();  _Vn1 = Vn1.V(); _parity = Vn.parity(); }
        void to_Lucas(arithmetic::LucasV& Vn, arithmetic::LucasV& Vn1) { Vn.V() = _Vn; Vn1.V() = _Vn1; Vn.arithmetic().init(Vn.V(), _parity, Vn); Vn1.arithmetic().init(Vn1.V(), !_parity, Vn1); }
        arithmetic::Giant& Vn() { return _Vn; }
        arithmetic::Giant& Vn1() { return _Vn1; }
        bool parity() { return _parity; }
        bool read(Reader& reader) override { int parity = 0; bool res = TaskState::read(reader) && reader.read(_Vn) && reader.read(_Vn1) && reader.read(parity); _parity = parity == 1; return res; }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_Vn); writer.write(_Vn1); writer.write(_parity ? 1 : 0); }

    private:
        arithmetic::Giant _Vn;
        arithmetic::Giant _Vn1;
        bool _parity;
    };
    class StrongCheckState : public TaskState
    {
    public:
        static const char TYPE = 11;
        StrongCheckState() : TaskState(TYPE) { }
        void set(int iteration, int recovery, const arithmetic::LucasV& Vn, const arithmetic::LucasV& Vn1, const arithmetic::LucasUV& D) { TaskState::set(iteration); _recovery = recovery; _Vn = Vn.V(); _Vn1 = Vn1.V(); _Vparity = Vn.parity(); _U = D.U(); _V = D.V(); _parity = D.parity(); }
        int recovery() { return _recovery; }
        void to_Lucas(arithmetic::LucasV& Vn, arithmetic::LucasV& Vn1, arithmetic::LucasUV& D) { Vn.V() = _Vn; Vn1.V() = _Vn1; Vn.arithmetic().init(Vn.V(), _Vparity, Vn); Vn1.arithmetic().init(Vn1.V(), !_Vparity, Vn1); D.U() = _U; D.V() = _V; D.arithmetic().init(D.U(), D.V(), _parity, true, D); }
        bool read(Reader& reader) override { int Vparity = 0, parity = 0; bool res = TaskState::read(reader) && reader.read(_recovery) && reader.read(_Vn) && reader.read(_Vn1) && reader.read(Vparity) && reader.read(_U) && reader.read(_V) && reader.read(parity); _Vparity = Vparity == 1; _parity = parity == 1; return res; }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_recovery); writer.write(_Vn); writer.write(_Vn1); writer.write(_Vparity ? 1 : 0); writer.write(_U); writer.write(_V); writer.write(_parity ? 1 : 0); }

    private:
        int _recovery;
        arithmetic::SerializedGWNum _Vn;
        arithmetic::SerializedGWNum _Vn1;
        bool _Vparity;
        arithmetic::SerializedGWNum _U;
        arithmetic::SerializedGWNum _V;
        bool _parity;
    };

public:
    template<class T>
    LucasUVMulFast(T&& exp, int count) : LucasVMul()
    {
        _exp = std::forward<T>(exp);
        Gerbicz_params((_exp.bitlen() + count - 1)/count, _L, _L2);
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging, int P, bool negativeQ)
    {
        _P = P;
        _negativeQ = negativeQ;
        init(input, gwstate, file, file_recovery, logging);
    }
    void init_state(State* state);

    State* state() { return _state_recovery.get(); }
    StrongCheckState* state_check() { return dynamic_cast<StrongCheckState*>(Task::state()); }
    arithmetic::Giant* result() override { if (state() == nullptr || state()->iteration() != iterations()) return nullptr; return &state()->Vn(); }

    void Gerbicz_params(int iters, int& L, int &L2);
    double cost() override { return _exp.bitlen()*2; }

    arithmetic::Giant& exp() { return _exp; }

protected:
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, File* file_recovery, Logging* logging);
    void write_state() override;
    void setup() override;
    void release() override;
    void execute() override;

    arithmetic::LucasUVArithmetic& arithmetic() { return *_arithmetic; }
    arithmetic::LucasUV& X() { return *_X; }
    arithmetic::LucasUV& R() { return *_R; }
    arithmetic::LucasUV& D() { return *_D; }

protected:
    arithmetic::Giant _exp;
    int _L;
    int _L2;

    int _P;
    bool _negativeQ;
    int _W;
    arithmetic::Giant _result;

    File* _file_recovery = nullptr;
    std::unique_ptr<State> _state_recovery;
    std::unique_ptr<State> _tmp_state_recovery;
    int _recovery_op = 0;

    std::unique_ptr<arithmetic::LucasUVArithmetic> _arithmetic;
    std::unique_ptr<arithmetic::LucasUV> _X;
    std::unique_ptr<arithmetic::LucasUV> _R;
    std::unique_ptr<arithmetic::LucasUV> _D;
};
