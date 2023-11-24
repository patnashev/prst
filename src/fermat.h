#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "exp.h"
#include "options.h"
#include "proof.h"

class Fermat
{
public:
    static const int AUTO = 0;
    static const int FERMAT = 1;
    static const int PROTH = 2;
    static const int POCKLINGTON = 3;

public:
    Fermat(int type, InputNum& input, Options& options, Logging& logging, Proof* proof);
    virtual ~Fermat() { }

    virtual void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof);

    int type() { return _type; }
    int a() { return _a; }
    bool success() { return _success; }
    bool prime() { return _prime; }
    std::string& res64() { return _res64; }
    arithmetic::Giant& result() { return *_task->result(); }

    CarefulExp* task_tail_simple() { return _task_tail_simple.get(); }
    CarefulExp* task_ak_simple() { return _task_ak_simple.get(); }
    CarefulExp* task_fermat_simple() { return _task_fermat_simple.get(); }
    MultipointExp* task() { return _task.get(); }

protected:
    int _type;
    int _a;
    Proof* _proof;

    std::unique_ptr<CarefulExp> _task_tail_simple;
    std::unique_ptr<CarefulExp> _task_ak_simple;
    std::unique_ptr<CarefulExp> _task_fermat_simple;
    std::unique_ptr<MultipointExp> _task;

    bool _success = false;
    bool _prime = false;
    std::string _res64;
};
