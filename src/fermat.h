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
    Fermat(InputNum& input, Params& params, Logging& logging, Proof* proof);

    void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof);

    bool Proth() { return _Proth; }
    int a() { return _a; }
    bool success() { return _success; }
    std::string& res64() { return _res64; }

    BaseExp* task_ak_simple() { return _task_ak_simple.get(); }
    BaseExp* task_ak() { return _task_ak.get(); }
    BaseExp* task() { return _task.get(); }

protected:

protected:
    bool _Proth = false;
    std::unique_ptr<InputNum> _input_k;
    std::unique_ptr<InputNum> _input_base2;
    int _a;

    std::unique_ptr<BaseExp> _task_ak_simple;
    std::unique_ptr<BaseExp> _task_ak;
    std::unique_ptr<BaseExp> _task;
    bool _success = false;
    std::string _res64;
};
