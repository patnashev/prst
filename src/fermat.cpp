
#include "gwnum.h"
#include "cpuid.h"
#include "fermat.h"
#include "integer.h"

using namespace arithmetic;

// LLR code for compatibility
uint32_t twopownmodm(uint32_t n, uint32_t m, uint32_t *order, uint32_t *nmodorder) {
    uint32_t tpnmodm, tp, i, work, mask = 1<<31;
    unsigned __int64 ltp;
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

Fermat::Fermat(InputNum& input, Params& params, Logging& logging, Proof* proof)
{
    if (input.b() == 2 && input.c() == 1)
        _Proth = true;
    else if (input.is_base2() && input.c() == 1 && proof == nullptr)
    {
        _input_k.reset(new InputNum());
        _input_base2.reset(new InputNum());
        input.to_base2(*_input_k, *_input_base2);
        _Proth = true;
        logging.info("%s = %s\n", input.display_text().data(), _input_base2->display_text().data());
    }

    _a = params.FermatBase ? params.FermatBase.value() : 3;
    if (Proth())
        _a = genProthBase(_input_base2 ? _input_base2->gk() : input.gk(), input.n());

    bool CheckGerbicz = true;// params.CheckGerbicz ? params.CheckGerbicz.value() : input.b() == 2;

    if (proof == nullptr && !CheckGerbicz)
    {
        FastExp* task;
        Giant exp;
        if (_input_base2)
            exp = _input_base2->gk() << (_input_base2->n() - 1);
        else
        {
            if (input.b() == 2)
                exp = input.gk() << (input.n() - (Proth() ? 1 : 0));
            else
                exp = input.gk()*power(input.gb(), input.n());
            exp += input.c() - 1;
        }
        _task.reset(task = new FastExp(std::move(exp)));
        logging.progress().add_stage(task->exp().bitlen());
        params.maxmulbyconst = _a;
    }
    else if (proof == nullptr)
    {
        Giant& b = _input_base2 ? _input_base2->gb() : input.gb();
        int n = (_input_base2 ? _input_base2->n() : input.n()) - (Proth() ? 1 : 0);
        int checks = params.GerbiczCount ? params.GerbiczCount.value() : 16;
        GerbiczCheckExp* task;
        _task.reset(task = params.GerbiczL ? new GerbiczCheckExp(b, n, checks, params.GerbiczL.value()) : new GerbiczCheckExp(b, n, checks));
        if (params.SlidingWindow)
            task->_W = params.SlidingWindow.value();

        if (!_input_k && input.k() != 1)
            _task_ak_simple.reset(new SlowExp(input.gk()));
        else if (_input_k)
        {
            if (_input_k->k() != 1)
                _task_ak_simple.reset(new SlowExp(_input_k->gk()));
            GerbiczCheckExp* task_ak;
            _task_ak.reset(task_ak = params.GerbiczL ? new GerbiczCheckExp(_input_k->gb(), _input_k->n(), checks, params.GerbiczL.value()) : new GerbiczCheckExp(_input_k->gb(), _input_k->n(), checks));
            if (params.SlidingWindow)
                task_ak->_W = params.SlidingWindow.value();
            logging.progress().add_stage(task_ak->cost());
        }

        logging.progress().add_stage(task->cost());
    }
    else
    {
        int n = input.n() - (Proth() ? 1 : 0);
        proof->calc_points(n, input, params, logging);
        CheckGerbicz = CheckGerbicz || params.ProofChecksPerPoint || params.ProofPointsPerCheck;
        MultipointExp* task;
        _task.reset(task = CheckGerbicz ? new GerbiczCheckMultipointExp(input.gb(), proof->points(), params.ProofPointsPerCheck ? params.ProofPointsPerCheck.value() : 1, params.ProofChecksPerPoint ? params.ProofChecksPerPoint.value() : 1, std::bind(&Proof::on_point, proof, std::placeholders::_1, std::placeholders::_2)) : new MultipointExp(input.gb(), proof->points(), std::bind(&Proof::on_point, proof, std::placeholders::_1, std::placeholders::_2)));
        GerbiczCheckMultipointExp* taskCheck = dynamic_cast<GerbiczCheckMultipointExp*>(task);
        if (CheckGerbicz && params.GerbiczL)
        {
            taskCheck->_L = params.GerbiczL.value();
            if (params.GerbiczL2)
                taskCheck->_L2 = params.GerbiczL2.value();
            else
            {
                taskCheck->_L2 = proof->M()*(params.ProofPointsPerCheck ? params.ProofPointsPerCheck.value() : 1)/(params.ProofChecksPerPoint ? params.ProofChecksPerPoint.value() : 1);
                taskCheck->_L2 -= taskCheck->_L2%taskCheck->_L;
            }
        }
        if (params.SlidingWindow)
            task->_W = params.SlidingWindow.value();
        logging.progress().add_stage(task->cost());
        logging.progress().add_stage(proof->cost());

        if (input.k() != 1)
            _task_ak_simple.reset(new SlowExp(input.gk()));
    }

    _task->set_error_check(!params.CheckNear || params.CheckNear.value(), params.Check && params.Check.value());
    if (_task_ak)
        _task_ak->set_error_check(!params.CheckNear || params.CheckNear.value(), params.Check && params.Check.value());
    if (_task_ak_simple)
        _task_ak_simple->set_error_check(false, true);
}

void Fermat::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof)
{
    File* ak_checkpoint = nullptr;
    File* ak_recoverypoint = nullptr;

    if (Proth())
        logging.info("Proth test of %s, a = %d.\n", (_input_base2 ? *_input_base2 : input).display_text().data(), _a);
    else
        logging.info("Fermat probabilistic test of %s, a = %d.\n", input.display_text().data(), _a);
    logging.report_param("a", _a);

    FastExp* taskFast = dynamic_cast<FastExp*>(_task.get());
    if (taskFast != nullptr)
    {
        taskFast->init(_input_base2 ? _input_base2.get() : &input, &gwstate, &file_checkpoint, &logging, _a);
        taskFast->run();
    }

    MultipointExp* taskMultipoint = dynamic_cast<MultipointExp*>(_task.get());
    if (taskMultipoint != nullptr)
    {
        GerbiczCheckMultipointExp* taskGerbiczCheck = dynamic_cast<GerbiczCheckMultipointExp*>(taskMultipoint);
        if (taskGerbiczCheck != nullptr)
            taskGerbiczCheck->init(_input_base2 ? _input_base2.get() : &input, &gwstate, &file_checkpoint, &file_recoverypoint, &logging);
        else
            taskMultipoint->init(_input_base2 ? _input_base2.get() : &input, &gwstate, &file_checkpoint, &logging);
        if (proof != nullptr)
            proof->init_state(taskMultipoint, gwstate, input, logging, _a);
        if (taskMultipoint->state() == nullptr)
        {
            Giant tmp;
            tmp = _a;
            if (_task_ak)
            {
                logging.set_prefix("");
                logging.info("Calculating a^(%s).\n", _input_k->display_text().data());
                GerbiczCheckMultipointExp* akGerbiczCheck = dynamic_cast<GerbiczCheckMultipointExp*>(_task_ak.get());
                if (akGerbiczCheck != nullptr)
                {
                    ak_checkpoint = file_checkpoint.add_child("ak", File::unique_fingerprint(gwstate.fingerprint, "ak"));
                    ak_recoverypoint = file_recoverypoint.add_child("ak", File::unique_fingerprint(gwstate.fingerprint, "ak"));
                    akGerbiczCheck->init(_input_k.get(), &gwstate, ak_checkpoint, ak_recoverypoint, &logging);
                    if (akGerbiczCheck->state() == nullptr)
                    {
                        if (_task_ak_simple)
                        {
                            static_cast<SlowExp*>(_task_ak_simple.get())->init(_input_k.get(), &gwstate, nullptr, &logging, tmp);
                            _task_ak_simple->run();
                            tmp = std::move(_task_ak_simple->state()->X());
                        }
                        akGerbiczCheck->init_state(new BaseExp::State(0, std::move(tmp)));
                    }
                    akGerbiczCheck->run();
                    logging.progress().next_stage();
                }
                tmp = std::move(_task_ak->state()->X());
            }
            else if (_task_ak_simple)
            {
                static_cast<SlowExp*>(_task_ak_simple.get())->init(&input, &gwstate, nullptr, &logging, tmp);
                _task_ak_simple->run();
                tmp = std::move(_task_ak_simple->state()->X());
            }

            if (proof != nullptr)
                proof->on_point(0, tmp);
            taskMultipoint->init_state(new BaseExp::State(0, std::move(tmp)));
        }
        else if (_input_k)
            logging.progress().skip_stage();

        taskMultipoint->run();
    }
    if (Proth())
    {
        _task->state()->X() += 1;
        _task->state()->X() %= *gwstate.N;
    }

    if (!Proth() && _task->state()->X() == 1)
    {
        _success = true;
        logging.result(_success, "%s is a probable prime. Time: %.1f s.\n", input.display_text().data(), _task->timer());
        logging.result_save(input.input_text() + " is a probable prime. Time: " + std::to_string((int)_task->timer()) + " s.\n");
    }
    else if (Proth() && _task->state()->X() == 0)
    {
        _success = true;
        logging.result(_success, "%s is prime! Time: %.1f s.\n", input.display_text().data(), _task->timer());
        logging.result_save(input.input_text() + " is prime! Time: " + std::to_string((int)_task->timer()) + " s.\n");
    }
    else
    {
        _res64 = _task->state()->X().to_res64();
        logging.result(_success, "%s is not prime. RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), _task->timer());
        logging.result_save(input.input_text() + " is not prime. RES64: " + _res64 + ", time: " + std::to_string((int)_task->timer()) + " s.\n");
    }

    logging.progress().next_stage();
    if (proof != nullptr)
        proof->run(input, gwstate, logging);

    file_checkpoint.clear();
    file_recoverypoint.clear();
    if (ak_checkpoint != nullptr)
        ak_checkpoint->clear();
    if (ak_recoverypoint != nullptr)
        ak_recoverypoint->clear();
}
