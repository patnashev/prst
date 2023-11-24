#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "options.h"
#include "lucasmul.h"

class Morrison
{
public:
    Morrison(InputNum& input, Options& options, Logging& logging);

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
        std::unique_ptr<LucasVMulFast> taskFactor;
        std::unique_ptr<LucasVMulFast> taskCheck;
    };

protected:
    bool _all_factors = false;
    bool _LLR = false;
    std::unique_ptr<LucasVMul> _task;
    std::unique_ptr<LucasVMulFast> _taskCheck;
    std::vector<FactorTask> _factor_tasks;
    int _P;
    bool _negQ;

    bool _success = false;
    bool _prime = false;
    std::string _res64;
};

