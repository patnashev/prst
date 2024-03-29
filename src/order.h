#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "exp.h"
#include "options.h"

class Order
{
public:
    Order(InputNum& a, InputNum& input, Options& options, Logging& logging);

    void run(InputNum& a, Options& options, InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging);

protected:
    void create_tasks(arithmetic::Giant& a, Options& options, Logging& logging, bool restart);
    bool on_point(int index, BaseExp::State* state);

    class FactorTask
    {
    public:
        FactorTask(arithmetic::Giant b_, int ord_, int n_) : b(b_), ord(ord_), n(n_) { }
        arithmetic::Giant b;
        int ord;
        int n;
        std::unique_ptr<CarefulExp> task_sub;
        std::unique_ptr<CarefulExp> task_factor;
    };

protected:
    int _sub = 30;
    std::vector<std::pair<arithmetic::Giant, int>> _factors;
    std::vector<std::pair<arithmetic::Giant, int>> _order;

    std::unique_ptr<MultipointExp> _task;
    std::string _task_exp_str;
    std::vector<std::unique_ptr<MultipointExp>> _tasks_smooth;
    std::unique_ptr<CarefulExp> _task_check;
    std::vector<FactorTask> _tasks;
    int _task_break;
};
