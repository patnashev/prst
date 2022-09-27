#pragma once

#include <list>
#include "fermat.h"

class Pocklington : public Fermat
{
public:
    Pocklington(InputNum& input, Params& params, Logging& logging, Proof* proof);

    void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof) override;

protected:
    std::list<std::pair<SlowExp,int>> _tasks;
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
