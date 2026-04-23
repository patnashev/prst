#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"

#include "prst.h"
#include "lucasmul.h"
#include "pocklington.h"

class Morrison : public Run
{
public:
    Morrison(InputNum& input, Options& options, Logging& logging);
    virtual ~Morrison() { }

    void run(arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging) override;

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
    std::unique_ptr<LucasVMul> _task;
    std::unique_ptr<LucasVMulFast> _taskCheck;
    std::vector<FactorTask> _factor_tasks;
    int _P;
    bool _negQ;
};

class MorrisonGeneric : public Run
{
public:
    MorrisonGeneric(InputNum& input, Options& options, Logging& logging);

    void run(arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging) override;

protected:
    std::unique_ptr<SubLogging> _logging;
    arithmetic::Giant _done;
    std::set<int> _done_factors;
    std::unique_ptr<FactorTree> _tree;
    std::vector<int> _dac_index;
    int _P;
    bool _negQ;
};
