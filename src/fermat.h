#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "exp.h"
#include "params.h"
#include "fastdc.h"

class Fermat
{
public:
    Fermat(InputNum& input, Params& params, FastDC* fastdc, Logging& logging);

    void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging);

    bool Proth() { return _Proth; }
    int a() { return _a; }
    bool success() { return _success; }

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
};
