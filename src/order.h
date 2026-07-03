#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "exp.h"
#include "prst.h"

class Order : public Run
{
public:
    Order(InputNum& input, Options& options, Logging& logging);

    void run(arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging) override;
    std::string& result() { return _result; }

protected:
    Order(const char* name, InputNum& input, Options& options) : Run(name, input, options) { }
    MultipointExp* create_smooth_task(arithmetic::Giant& base, int power);
    void create_tasks(arithmetic::Giant& a, Logging& logging, bool restart);
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
    std::string _result;
};

class FermatDivisor : public Order
{
public:
    FermatDivisor(InputNum& input, Options& options, Logging& logging);

    void run(arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging) override;

protected:
    class Base
    {
    public:
        Base(int b)
        {
            base = b;
            str = std::to_string(b);
        }

        arithmetic::Giant base;
        std::string str;
        int power;
        std::unique_ptr<MultipointExp> task;
        arithmetic::Giant sub_val;
        arithmetic::Giant exp;
        arithmetic::Giant val;
    };

    void create_task(Logging& logging, Base& base);
    void run_task(arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Base& base);

protected:
    std::vector<Base> _bases;
};
