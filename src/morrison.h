#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "options.h"
#include "lucasmul.h"
#include "pocklington.h"

class Morrison
{
protected:
    Morrison() { }
public:
    Morrison(InputNum& input, Options& options, Logging& logging);

    virtual void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging);

    bool success() { return _success; }
    bool prime() { return _prime; }
    std::string& res64() { return _res64; }
    bool is_LLR() { return _LLR; }

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

class MorrisonGeneric : public Morrison
{
public:
    MorrisonGeneric(InputNum& input, Options& options, Logging& logging);

    virtual void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging) override;

protected:
    Options& _options;
    std::unique_ptr<SubLogging> _logging;

    arithmetic::Giant _done;
    std::set<int> _done_factors;
    std::unique_ptr<FactorTree> _tree;
    std::vector<int> _dac_index;
};
