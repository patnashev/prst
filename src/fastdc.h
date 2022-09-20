#pragma once

#include "arithmetic.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "params.h"

class FastDC
{
public:
    static const int NO_OP = 0;

public:
    FastDC(int op, InputNum& input, Params& params) : _op(op)
    {
    }

    int op() { return _op; }

protected:

protected:
    int _op;
};
