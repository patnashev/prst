#pragma once

#include "fermat.h"

class Pocklington : public Fermat
{
public:
    Pocklington(InputNum& input, Options& options, Logging& logging, Proof* proof);

    void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof) override;

protected:
    class FactorTask
    {
    public:
        FactorTask(int i) : index(i) { }
        int index;
        std::unique_ptr<CarefulExp> taskFactor;
        std::unique_ptr<CarefulExp> taskCheck;
    };

protected:
    bool _all_factors = false;
    std::vector<FactorTask> _tasks;
    arithmetic::Giant _done;
};
