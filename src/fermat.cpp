
#include <cmath>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "fermat.h"
#include "integer.h"

using namespace arithmetic;

// LLR code for compatibility
uint32_t twopownmodm(uint32_t n, uint32_t m, uint32_t *order, uint32_t *nmodorder) {
    uint32_t tpnmodm, tp, i, work, mask = 1<<31;
    uint64_t ltp;
    tpnmodm = 0;				// This function computes 2^n modulo m
    if (!(m&1)||(m==1)) {		// It returns this value, and also :
        *order = 0;				// The order of 2 modulo m, and the remainder of n modulo this order.
        *nmodorder = n;
    }
    if (m==1)
        return (tpnmodm);
    tp = 1;
    for (i = 1; i<m; i++) {
        tp <<= 1;
        if (tp >= m) tp -= m;	// Modular reduction
        if (i==n) tpnmodm = tp;	// If m is even or n < order of 2 modulo m the calculus is completed here.
        if (tp==1) {			// But continue to compute the order of 2
            *order = i;
            break;
        }
    }
    if (*order)
        *nmodorder = n%*order;	// Compute the remainder
    if (!tpnmodm) {				// If n >= order, continue
        work = *nmodorder;
        if (!work)				// n is a multiple of the order...
            return (1);
        while (!(work&mask))
            mask >>= 1;
        ltp = 1;				// init the result
        while (mask) {
            ltp *= ltp;			// square the result
            ltp %= m;			// modular reduction
            if (work&mask) {	// test the current bit of the exponent
                ltp <<= 1;		// multiply the result by the base
                if (ltp >= m) ltp -= m;
            }
            mask >>= 1;			// shift the mask
        }

        /*	    for (i=1, tp=1; i<=*nmodorder; i++) {
                    tp <<= 1;
                    if (tp >= m) tp -= m;
                } */
        tpnmodm = (uint32_t)ltp;	// 2^n modulo m == 2^(n%order) modulo m
    }
    return tpnmodm;
}

int genProthBase(Giant& k, uint32_t n) {
    uint32_t Nmodp, kmodp, p, tpnmp, orderp, nmodorderp, kw;
    int jNp;

    //	Return the least prime p such Jacobi (N, p) = -1

    if (k.size() == 1 && n < 3) {	//	Eliminate some trivial cases
        kw = k.data()[0];
        if (n == 1 && kw == 1)
            return (2);
        else if (n == 2)
            return (2);
        else
            return (-1);
    }
    else {							// General case
        for (p = 3; p<=2147483647; p += 2) {
            if (!is_prime(p))
                continue;
            kmodp = k%p;
            if (!kmodp)
                continue;
            tpnmp = twopownmodm(n, p, &orderp, &nmodorderp);
            Nmodp = (kmodp*tpnmp+1)%p;
            if (!Nmodp) {
                return (-(int)p);
            }
            if ((jNp = jacobi(Nmodp, p)) > 1) {
                return (-jNp);
            }
            if (jNp != -1)
                continue;
            return (p);
        }
        return (-1);
    }
}
// END LLR code

Fermat::Fermat(int type, InputNum& input, Params& params, Logging& logging, Proof* proof)
{
    if ((type == AUTO || type == PROTH) && input.b() == 2 && input.c() == 1)
        _type = PROTH;
    else if (type == AUTO)
        _type = FERMAT;
    else
        _type = type;

    _a = params.FermatBase ? params.FermatBase.value() : 3;
    if (_type == PROTH)
    {
        _a = genProthBase(input.gk(), input.n());
        if (_a < 0)
        {
            logging.result(false, "%s is not prime, divisible by %d.\n", input.display_text().data(), -_a);
            logging.result_save(input.input_text() + " is not prime, divisible by " + std::to_string(-_a) + ".\n");
        }
    }
    _n = input.n() - (_type == PROTH || _type == POCKLINGTON ? 1 : 0);

    bool CheckStrong = params.CheckStrong ? params.CheckStrong.value() : false;
    auto on_point = std::bind(&Fermat::on_point, this, std::placeholders::_1, std::placeholders::_2);
    Giant exp;
    if ((proof == nullptr && (!CheckStrong || input.b() != 2)) || (proof != nullptr && proof->Li()))
    {
        if (input.b() == 2)
            exp = input.gk() << _n;
        else
            exp = input.gk()*power(input.gb(), _n);
        exp += input.c() - 1;
    }

    if (proof == nullptr && !CheckStrong)
    {
        _task.reset(new FastExp(std::move(exp)));
        logging.progress().add_stage(_task->exp().bitlen());
        params.maxmulbyconst = _a;
        if (_type == PROTH || _type == POCKLINGTON)
            _task_b_simple.reset(new CarefulExp(input.gb()));
    }
    else if (proof == nullptr)
    {
        int checks = params.StrongCount ? params.StrongCount.value() : 16;
        if (input.b() == 2)
        {
            _task.reset(new GerbiczCheckExp(input.gb(), _n, checks, on_point, params.StrongL ? params.StrongL.value() : 0));
            if (_type == PROTH || _type == POCKLINGTON)
                _task->points().push_back(_n + 1);
            _points = _task->points();

            if (input.c() != 1)
                _task_tail_simple.reset(new CarefulExp(abs(input.c() - 1)));
            if (input.k() != 1)
                _task_ak_simple.reset(new CarefulExp(input.gk()));
        }
        else
        {
            _task.reset(new FastLiCheckExp(std::move(exp), checks, params.StrongL ? params.StrongL.value() : 0));
            params.maxmulbyconst = _a;
            if (_type == PROTH || _type == POCKLINGTON)
                _task_b_simple.reset(new CarefulExp(input.gb()));
        }
        logging.progress().add_stage(_task->cost());
    }
    else
    {
        if (proof->Li())
            proof->calc_points(exp.bitlen() - 1, input, params, logging);
        else
            proof->calc_points(_n, input, params, logging);
        int L = 0, L2 = 0;
        if (CheckStrong)
        {
            L2 = proof->M()*(params.ProofPointsPerCheck ? params.ProofPointsPerCheck.value() : 1)/(params.ProofChecksPerPoint ? params.ProofChecksPerPoint.value() : 1);
            if (params.StrongL)
            {
                L = params.StrongL.value();
                if (params.StrongL2)
                    L2 = params.StrongL2.value();
                else
                    L2 -= L2%L;
            }
            else
                StrongCheckMultipointExp::Gerbicz_params(L2, 1.0, L, L2);
        }

        if (proof->Li())
        {
            auto on_point = std::bind(&Proof::on_point, proof, std::placeholders::_1, std::placeholders::_2);
            if (!CheckStrong)
                _task.reset(new MultipointExp(std::move(exp), false, proof->points(), on_point));
            else
                _task.reset(new StrongCheckMultipointExp(std::move(exp), false, proof->points(), L, L2, on_point));
            params.maxmulbyconst = _a;
            if (_type == PROTH || _type == POCKLINGTON)
                _task_b_simple.reset(new CarefulExp(input.gb()));
        }
        else
        {
            _points = proof->points();
            if (_points.back() != _n)
                _points.push_back(_n);
            if (_type == PROTH || _type == POCKLINGTON)
                _points.push_back(_n + 1);

            if (!CheckStrong)
                _task.reset(new MultipointExp(input.gb(), true, _points, on_point));
            else
                _task.reset(new StrongCheckMultipointExp(input.gb(), true, _points, L, L2, on_point));
            if (params.SlidingWindow)
                _task->_W = params.SlidingWindow.value();

            if (input.c() != 1)
                _task_tail_simple.reset(new CarefulExp(abs(input.c() - 1)));
            if (input.k() != 1)
                _task_ak_simple.reset(new CarefulExp(input.gk()));
        }
        logging.progress().add_stage(_task->cost());
        logging.progress().add_stage(proof->cost());
    }

    _task->set_error_check(!params.CheckNear || params.CheckNear.value(), params.Check && params.Check.value());
    if (_task_tail_simple)
        _task_tail_simple->set_error_check(false, true);
    if (_task_ak_simple)
        _task_ak_simple->set_error_check(false, true);
    if (_task_b_simple)
        _task_b_simple->set_error_check(false, true);
}

bool Fermat::on_point(int index, arithmetic::Giant& X)
{
    bool ret = false;
    if (_proof != nullptr && index <= _proof->count())
        ret = _proof->on_point(index, X);
    if ((_type == PROTH || _type == POCKLINGTON) && _points[index] == _n)
        _Xm1 = X;
    return ret;
}

void Fermat::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof)
{
    _proof = proof;
    _success = false;
    _res64 = "";
    _Xm1 = Giant();
    Giant ak;
    ak = _a;
    if (_a < 0)
        return;

    if (type() == PROTH)
        logging.info("Proth test of %s, a = %d, complexity = %d.\n", input.display_text().data(), _a, (int)logging.progress().cost_total());
    else if (type() == FERMAT)
        logging.info("Fermat probabilistic test of %s, a = %d, complexity = %d.\n", input.display_text().data(), _a, (int)logging.progress().cost_total());
    logging.report_param("a", _a);
    if (gwstate.information_only)
        exit(0);

    Giant tail;
    if (_task_tail_simple)
    {
        if (input.c() == -1 && _a < 46341)
            tail = _a*_a;
        else
        {
            _task_tail_simple->init(&input, &gwstate, &logging, ak);
            _task_tail_simple->run();
            tail = std::move(_task_tail_simple->state()->X());
        }
        if (input.c() < 0)
        {
            tail.inv(*gwstate.N);
            if (input.c() == -1 && _a < 46341)
                GWASSERT(tail*(_a*_a)%*gwstate.N == 1);
        }
    }

    StrongCheckMultipointExp* taskCheck = dynamic_cast<StrongCheckMultipointExp*>(_task.get());
    if (taskCheck != nullptr && taskCheck->smooth())
        taskCheck->init_smooth(&input, &gwstate, &file_checkpoint, &file_recoverypoint, &logging, std::move(tail));
    else if (taskCheck != nullptr)
        taskCheck->init_small(&input, &gwstate, &file_checkpoint, &file_recoverypoint, &logging, _a, std::move(tail));
    else if (_task->smooth())
        _task->init_smooth(&input, &gwstate, &file_checkpoint, &logging, std::move(tail));
    else
        _task->init_small(&input, &gwstate, &file_checkpoint, &logging, _a, std::move(tail));
    if (proof != nullptr)
        proof->init_state(_task.get(), gwstate, input, logging, _a);
    if (_task->smooth() && (_task->state() == nullptr || _task->state()->iteration() > _n))
    {
        if (_task_ak_simple)
        {
            _task_ak_simple->init(&input, &gwstate, &logging, std::move(ak));
            _task_ak_simple->run();
            ak = std::move(_task_ak_simple->state()->X());
        }

        if (proof != nullptr)
        {
            logging.progress().update(0, (int)gwstate.handle.fft_count/2);
            logging.progress_save();
            proof->on_point(0, ak);
        }
        _task->init_state(new BaseExp::State(0, std::move(ak)));
        if (proof != nullptr)
            _task->state()->set_written();
    }
    if (_task->smooth() && (_type == PROTH || _type == POCKLINGTON) && _task->state()->iteration() == _n)
        _Xm1 = _task->state()->X();

    _task->run();

    if (_task_b_simple)
    {
        _Xm1 = std::move(_task->state()->X());
        _task_b_simple->init(&input, &gwstate, &logging, _Xm1);
        _task_b_simple->run();
        _task->state()->X() = std::move(_task_b_simple->state()->X());
    }

    if (type() == PROTH)
        _Xm1 += 1;

    if (type() == PROTH && (_Xm1 == 0 || _Xm1 == *gwstate.N))
    {
        _success = true;
        _prime = true;
        logging.result(_prime, "%s is prime! Time: %.1f s.\n", input.display_text().data(), logging.progress().time_total());
        logging.result_save(input.input_text() + " is prime! Time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
    }
    else if (type() == PROTH || _task->state()->X() != 1)
    {
        if (type() == PROTH)
        {
            _res64 = _Xm1.to_res64();
            _Xm1 -= 1;
        }
        else
            _res64 = _task->state()->X().to_res64();
        logging.result(_prime, "%s is not prime. RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), logging.progress().time_total());
        logging.result_save(input.input_text() + " is not prime. RES64: " + _res64 + ", time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
    }
    if (!_prime && _task->state()->X() == 1)
    {
        _success = true;
        logging.result(type() != PROTH && type() != POCKLINGTON, "%s is a probable prime. Time: %.1f s.\n", input.display_text().data(), logging.progress().time_total());
        logging.result_save(input.input_text() + " is a probable prime. Time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
    }

    logging.progress().next_stage();
    if (proof != nullptr)
        proof->run(input, gwstate, logging, _success ? nullptr : &result());

    file_checkpoint.clear();
    file_recoverypoint.clear();
}
