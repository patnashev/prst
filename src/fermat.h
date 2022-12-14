#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "exp.h"
#include "params.h"
#include "proof.h"

class Fermat
{
public:
    static const int AUTO = 0;
    static const int FERMAT = 1;
    static const int PROTH = 2;
    static const int POCKLINGTON = 3;

public:
    Fermat(int type, InputNum& input, Params& params, Logging& logging, Proof* proof);
    virtual ~Fermat() { }

    virtual void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof);

    int type() { return _type; }
    int a() { return _a; }
    bool success() { return _success; }
    std::string& res64() { return _res64; }
    arithmetic::Giant& result() { return _Xm1.empty() ? _task->state()->X() : _Xm1; }

    CarefulExp* task_tail_simple() { return _task_tail_simple.get(); }
    CarefulExp* task_ak_simple() { return _task_ak_simple.get(); }
    CarefulExp* task_b_simple() { return _task_b_simple.get(); }
    MultipointExp* task() { return _task.get(); }

protected:
    bool on_point(int index, arithmetic::Giant& X);

protected:
    int _type;
    int _a;
    int _n;
    Proof* _proof;
    std::vector<int> _points;
    arithmetic::Giant _Xm1;

    std::unique_ptr<CarefulExp> _task_tail_simple;
    std::unique_ptr<CarefulExp> _task_ak_simple;
    std::unique_ptr<CarefulExp> _task_b_simple;
    std::unique_ptr<MultipointExp> _task;

    bool _success = false;
    std::string _res64;
};
