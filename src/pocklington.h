#pragma once

#include "fermat.h"

class Pocklington : public Fermat
{
public:
    Pocklington(InputNum& input, Params& params, Logging& logging, Proof* proof);

    void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof) override;

protected:
    void on_finish(InputNum& input, arithmetic::GWState& gwstate, Logging& logging) override;

protected:
    std::list<SlowExp> _tasks;
    std::unique_ptr<SlowExp> _task_fermat_simple;
    arithmetic::Giant _Xm1;
};

template<class IT>
class Product : public BaseExp
{
public:
    Product(IT first, IT last) : BaseExp(), _first(first), _last(last)
    {
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging);

protected:
    void execute() override;

protected:
    IT _first;
    IT _last;
};
