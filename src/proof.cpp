
#include <cmath>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "proof.h"
#include "md5.h"
#include "integer.h"

using namespace arithmetic;

Proof::Proof(int op, int count, InputNum& input, Params& params, File& file_cert, Logging& logging, std::optional<bool> forceLi) : _op(op), _count(count)
{
    if ((op == SAVE || op == BUILD) && (count & (count - 1)) != 0)
    {
        logging.error("proof count is not a power of 2.\n");
        exit(0);
    }

    bool CheckStrong = params.CheckStrong ? params.CheckStrong.value() : false;
    _Li = forceLi ? forceLi.value() : input.b() != 2;

    if (op == SAVE)
        _task.reset(new ProofSave(*this));
    if (op == BUILD)
        _task.reset(new ProofBuild(*this, params.ProofSecuritySeed));
    if (op == CERT)
    {
        Certificate cert;
        if (!file_cert.read(cert) || (Li() && cert.a_power() == 0))
        {
            logging.error("Invalid certificate file.\n");
            exit(0);
        }
        file_cert.free_buffer();
        _M = cert.power();
        _r_count = std::move(cert.X());
        _r_0 = std::move(cert.a_base());
        int checks = params.StrongCount ? params.StrongCount.value() : 1;

        Giant b = input.gb();
        if (Li())
            b = 2;
        MultipointExp* task;
        if (!CheckStrong)
            _task.reset(task = new SmoothExp(b, _M));
        else
            _task.reset(task = new GerbiczCheckExp(b, _M, checks, nullptr, params.StrongL ? params.StrongL.value() : 0));
        if (params.SlidingWindow)
            task->_W = params.SlidingWindow.value();

        if (Li())
        {
            if (!CheckStrong)
                _taskA.reset(task = new SlidingWindowExp(std::move(cert.a_power())));
            else
                _taskA.reset(task = new LiCheckExp(std::move(cert.a_power()), checks, params.StrongL ? params.StrongL.value() : 0));
            if (params.SlidingWindow)
                task->_W = params.SlidingWindow.value();
            logging.progress().add_stage(_taskA->cost());
        }
        logging.progress().add_stage(_task->cost());
    }
    if ((op == BUILD && (!params.RootOfUnityCheck || params.RootOfUnityCheck.value())) || op == ROOT)
    {
        Giant exp;
        if (input.c() == 1)
        {
            exp = input.gk();
            int security = params.RootOfUnitySecurity ? params.RootOfUnitySecurity.value() : 64;
            for (auto it = input.b_factors().begin(); it != input.b_factors().end(); it++)
                if (it->first == 2)
                    exp <<= security;
                else
                    for (double bitlen = 0, log2factor = log2(it->first); bitlen < security; bitlen += log2factor)
                        exp *= it->first;
        }
        else
        {
            exp = 1;
            int security = params.RootOfUnitySecurity ? params.RootOfUnitySecurity.value() : 24;
            logging.info("Factorizing N-1 for roots of unity check");
            double timer = getHighResTimer();
            std::vector<int> factors = input.factorize_minus1(security);
            timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
            logging.info(", time: %.1f s.\n", timer);
            for (auto it = factors.begin(); it != factors.end(); it++)
                exp *= *it;
        }
        _taskRoot.reset(new CarefulExp(std::move(exp)));
    }

    if (_task)
        _task->set_error_check(params.CheckNear && params.CheckNear.value(), !params.Check || params.Check.value());
    if (_taskA)
        _taskA->set_error_check(params.CheckNear && params.CheckNear.value(), !params.Check || params.Check.value());
    if (_taskRoot)
        _taskRoot->set_error_check(params.CheckNear && params.CheckNear.value(), !params.Check || params.Check.value());
}

void Proof::calc_points(int iterations, InputNum& input, Params& params, Logging& logging)
{
    if (params.StrongCount)
        if (params.StrongCount.value() > _count)
            params.ProofChecksPerPoint = params.StrongCount.value()/_count;
        else
            params.ProofPointsPerCheck = _count/params.StrongCount.value();
    _points_per_check = 1;
    if ((input.b() != 2 || input.c() != 1) && params.ProofPointsPerCheck && !Li())
    {
        int i;
        _points_per_check = params.ProofPointsPerCheck.value();
        int iters = iterations*_points_per_check/_count;
        if (!params.StrongL)
        {
            int L = _points_per_check*(int)std::sqrt(iters);
            int L2 = iters - iters%L;
            for (i = L + _points_per_check; i*i < 2*iters*_points_per_check*_points_per_check; i += _points_per_check)
                if (L2 < iters - iters%i)
                {
                    L = i;
                    L2 = iters - iters%i;
                }
            params.StrongL = L/_points_per_check;
            _M = L2/_points_per_check;
        }
        else
            _M = (iters - iters%(params.StrongL.value()*_points_per_check))/_points_per_check;
        _points.reserve(_count + 2);
        for (i = 0; i <= _count; i++)
            _points.push_back(i*_M);
    }
    else
    {
        _points.reserve(_count + 1);
        _points.push_back(0);
        for (int i = 1; i < _count; i++)
        {
            _M = iterations;
            _points.push_back(0);
            for (int j = _count/2; j > 0 && (i & (j*2 - 1)) != 0; j >>= 1)
            {
                _M /= 2;
                if ((i & j) != 0)
                    _points.back() += _M;
                if ((iterations & (_count/j/2)) != 0)
                    _points.back()++;
            }
        }
        _points.push_back(iterations);
    }
}

void Proof::init_files(File* file_point, File* file_product, File* file_cert)
{
    int i;
    _file_points.clear();
    _file_points.reserve(_count + 1);
    for (i = 0; i <= _count; i++)
        _file_points.push_back(file_point->add_child(std::to_string(i), file_point->fingerprint()));
    _file_products.clear();
    if (file_product != nullptr)
        for (i = 0; (1 << i) < _count; i++)
            _file_products.push_back(file_product->add_child(std::to_string(i), file_product->fingerprint()));
    _file_cert = file_cert;
}

void Proof::init_state(MultipointExp* task, arithmetic::GWState& gwstate, InputNum& input, Logging& logging, int a)
{
    std::unique_ptr<BaseExp::State> state(new BaseExp::State());

    if (op() == BUILD)
    {
        if (_points.size() == count() + 2)
            logging.info("Building certificate from %d products, %d iterations tail.\n", depth(), _points[count() + 1] - _points[count()]);
        else
            logging.info("Building certificate from %d products.\n", depth());
        logging.set_prefix(input.display_text() + " ");

        if (!Li())
        {
            read_point(0, *state, logging);
            double timer = getHighResTimer();
            Giant tmp;
            tmp = a;
            tmp.arithmetic().powermod(tmp, input.gk(), *gwstate.N, tmp);
            if (tmp != state->X())
            {
                logging.error("invalid a^k.\n");
                throw TaskAbortException();
            }
            timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
            logging.info("a^k is correct. Time: %.1f s.\n", timer);
            _r_0 = std::move(state->X());
        }
        else
        {
            _r_0 = a;
            _r_exp = &task->exp();
        }

        read_point(_count, *state, logging);
        _r_count = state->X();

        task->set_error_check(false, true);
        task->init_state(state.release());
        return;
    }

    logging.info("Saving %d proof points.\n", count());
    int point = _count;
    while (point >= 0)
    {
        point -= point%_points_per_check;
        if (_file_points[point]->read(*state) && state->iteration() == _points[point])
        {
            if (task->state() == nullptr || task->state()->iteration() < state->iteration() || (point < _count && task->state()->iteration() >= _points[point + 1]))
                task->init_state(state.release());
            return;
        }
        point--;
    }
    if (task->state() != nullptr)
        task->init_state(nullptr);
}

void Proof::read_point(int index, TaskState& state, Logging& logging)
{
    if (!_file_points[index]->read(state) || state.iteration() != _points[index])
    {
        logging.error("%s is missing or corrupt.\n", _file_points[index]->filename().data());
        throw TaskAbortException();
    }
    if (!_cache_points)
        _file_points[index]->free_buffer();
}

void Proof::read_product(int index, TaskState& state, Logging& logging)
{
    if (!_file_products[index]->read(state) || state.iteration() != index)
    {
        logging.error("%s is missing or corrupt.\n", _file_products[index]->filename().data());
        throw TaskAbortException();
    }
    _file_products[index]->free_buffer();
}

bool Proof::on_point(int index, arithmetic::Giant& X)
{
    if (index > _count)
        return false;
    BaseExp::State state(_points[index], std::move(X));
    _file_points[index]->write(state);
    if (!_cache_points)
        _file_points[index]->free_buffer();
    X = std::move(state.X());
    return true;
}

void Proof::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging)
{
    double timer = 0;
    Giant tail;
    if (Li())
    {
        File* file_checkpoint_a(file_checkpoint.add_child("a", file_checkpoint.fingerprint()));
        File* file_recoverypoint_a(file_recoverypoint.add_child("a", file_recoverypoint.fingerprint()));
        logging.info("Verifying certificate of %s, %d+%d iterations.\n", input.display_text().data(), _taskA->exp().bitlen() - 1, _M);

        MultipointExp* taskA = dynamic_cast<MultipointExp*>(_taskA.get());
        LiCheckExp* taskACheck = dynamic_cast<LiCheckExp*>(_taskA.get());
        if (taskACheck != nullptr)
            taskACheck->init(&input, &gwstate, file_checkpoint_a, file_recoverypoint_a, &logging, std::move(_r_0));
        else
            taskA->init(&input, &gwstate, file_checkpoint_a, &logging, std::move(_r_0));
        _taskA->run();
        timer = _taskA->timer();
        tail = std::move(_taskA->state()->X());

        logging.progress().next_stage();
        file_checkpoint_a->clear();
        file_recoverypoint_a->clear();
    }
    else
        logging.info("Verifying certificate of %s, %d iterations.\n", input.display_text().data(), _M);

    MultipointExp* task = dynamic_cast<MultipointExp*>(_task.get());
    if (task != nullptr)
    {
        GerbiczCheckExp* taskCheck = dynamic_cast<GerbiczCheckExp*>(_task.get());
        if (taskCheck != nullptr)
            taskCheck->init(&input, &gwstate, &file_checkpoint, &file_recoverypoint, &logging, std::move(tail));
        else
            task->init_smooth(&input, &gwstate, &file_checkpoint, &logging, std::move(tail));
        if (task->state() == nullptr)
            task->init_state(new BaseExp::State(0, std::move(_r_count)));
    }
    _task->run();
    timer += _task->timer();

    _res64 = _task->state()->X().to_res64();
    logging.result(false, "%s certificate RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), timer);
    logging.result_save(input.input_text() + " certificate RES64: " + _res64 + ", time: " + std::to_string((int)timer) + " s.\n");

    logging.progress().next_stage();
    file_checkpoint.clear();
    file_recoverypoint.clear();
}

void Proof::run(InputNum& input, arithmetic::GWState& gwstate, Logging& logging, Giant* X)
{
    CarefulExp* taskRoot = dynamic_cast<CarefulExp*>(_taskRoot.get());
    if (taskRoot != nullptr && X != nullptr)
    {
        taskRoot->init(&input, &gwstate, nullptr, &logging, std::move(*X));
        taskRoot->run();
        if (taskRoot->state()->X() == 1)
        {
            logging.error("%s roots of unity check failed.\n", input.display_text().data());
            throw TaskAbortException();
        }
    }
    ProofSave* taskSave = dynamic_cast<ProofSave*>(_task.get());
    if (taskSave != nullptr)
    {
        taskSave->init(&input, &gwstate, &logging);
        taskSave->run();
        logging.progress().next_stage();

        _res64 = taskSave->state()->X().to_res64();
        logging.info("%s compressed %d points to %d products, time: %.1f s.\n", input.display_text().data(), _count, _file_products.size(), _task->timer());
        logging.result(false, "%s raw certificate RES64: %s.\n", input.display_text().data(), _res64.data());
        logging.result_save(input.input_text() + " raw certificate RES64: " + _res64 + ", time: " + std::to_string((int)_task->timer()) + " s.\n");
    }
    ProofBuild* taskBuild = dynamic_cast<ProofBuild*>(_task.get());
    if (taskBuild != nullptr)
    {
        taskBuild->init(&input, &gwstate, &logging);
        taskBuild->run();
        logging.progress().next_stage();

        if (taskBuild->security())
        {
            _res64 = taskBuild->raw_res64();
            logging.result(false, "%s raw certificate RES64: %s.\n", input.display_text().data(), _res64.data());
            logging.result_save(input.input_text() + " raw certificate RES64: " + _res64 + ".\n");
        }
        _res64 = taskBuild->state()->X().to_res64();
        logging.result(false, "%s certificate RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), _task->timer());
        logging.result_save(input.input_text() + " certificate RES64: " + _res64 + ", time: " + std::to_string((int)_task->timer()) + " s.\n");
    }
}

double Proof::cost()
{
    if (op() == SAVE)
        return count()*depth()*16;
    if (op() == BUILD)
        return (depth() + (static_cast<ProofBuild*>(_task.get())->security() ? 1 : 0))*2*96;
    return 0;
}

void ProofSave::init(InputNum* input, arithmetic::GWState* gwstate, Logging* logging)
{
    BaseExp::init(input, gwstate, nullptr, nullptr, logging, _proof.count());
    _state_update_period = 1;
    _logging->set_prefix(input->display_text() + " ");
}

void ProofSave::read_point(int index, TaskState& state)
{
    _proof.read_point(index, state, *_logging);
}

void hash_giant(Giant& gin, Giant& gout)
{
    MD5_CTX context;
    MD5Init(&context);
    MD5Update(&context, (unsigned char *)gin.data(), gin.size()*4);
    MD5Final((unsigned char *)gout.data(), &context);
    gout.arithmetic().init(gout.data(), 2, gout);
}

void hash_giants(uint32_t fingerprint, Giant& gin1, Giant& gin2, Giant& gout)
{
    MD5_CTX context;
    MD5Init(&context);
    MD5Update(&context, (unsigned char *)&fingerprint, 4);
    MD5Update(&context, (unsigned char *)gin1.data(), gin1.size()*4);
    MD5Update(&context, (unsigned char *)gin2.data(), gin2.size()*4);
    MD5Final((unsigned char *)gout.data(), &context);
    gout.arithmetic().init(gout.data(), 2, gout);
}

void make_prime(Giant& g, int limit = 1000)
{
    g.data()[0] |= 1;
    if (g.size() < 2)
        g.data()[1] = 1;
    g.arithmetic().init(g.data(), 2, g);
    while (true)
    {
        PrimeIterator it = PrimeIterator::get();
        for (; *it < limit && g%(*it) != 0; it++);
        if (*it >= limit)
            break;
        g += 2;
    }
}

void exp_gw(GWArithmetic& gw, Giant& exp, GWNum& X, GWNum& X0, int options)
{
    int len = exp.bitlen() - 1;
    for (int bit = 0; bit < len; bit++)
    {
        gw.square(X, X, options);
        if (exp.bit(len - bit - 1))
            gw.mul(X0, X, X, options);
    }
}

void ProofSave::execute()
{
    int t, i, j, k;
    Proof::Product state_d;
    GWNum Y(gw());
    GWNum D(gw());
    GWNum T(gw());
    std::vector<GWNum> tree;
    std::vector<Giant> h;

    t = _proof.depth();
    tree.reserve(t);
    h.reserve(t);

    if (state() == nullptr)
    {
        _state.reset(new State());
        read_point(_proof.count(), *state());
        static_cast<TaskState*>(state())->set(0);
    }
    Y = state()->X();

    for (i = state()->iteration(); i < t; i++, commit_execute<State>(i, Y))
    {
        if (_proof.file_products()[i]->read(state_d))
        {
            _proof.file_products()[i]->free_buffer();
            D = state_d.X();
        }
        else
        {
            state_d.mimic_type(State::TYPE);
            if (i == 0)
            {
                read_point(_proof.count()/2, state_d);
                static_cast<TaskState&>(state_d).set(0);
                D = state_d.X();
            }
            else
            {
                while (tree.size() < i)
                    tree.emplace_back(gw());
                for (j = 0; j < (1 << i); j++)
                {
                    k = (1 + j*2) << (t - i - 1);
                    read_point(k, state_d);
                    D = state_d.X();

                    for (k = 1; k <= i; k++)
                    {
                        if (Task::abort_flag())
                            throw TaskAbortException();
                        if ((j & (1 << (k - 1))) == 0)
                        {
                            gw().fft(D, tree[i - k]);
                            break;
                        }
                        else
                        {
                            exp_gw(gw(), h[i - k], T = tree[i - k], tree[i - k], GWMUL_STARTNEXTFFT);
                            gw().mul(T, D, D, GWMUL_STARTNEXTFFT_IF(j + 1 != (1 << i) || k != i));
                        }
                    }
                }
                state_d.set(i, D);
            }

            state_d.mimic_type(Proof::Product::TYPE);
            _proof.file_products()[i]->write(state_d);
            _proof.file_products()[i]->free_buffer();
        }

        h.emplace_back(GiantsArithmetic::default_arithmetic(), 4);
        hash_giants(_gwstate->fingerprint, state()->X(), state_d.X(), h[i]);
        make_prime(h[i]);

        gw().fft(D, D);
        exp_gw(gw(), h[i], D, T = D, GWMUL_STARTNEXTFFT);
        gw().mul(D, Y, Y, 0);
    }
    tree.clear();

    done();
}

void ProofBuild::init(InputNum* input, arithmetic::GWState* gwstate, Logging* logging)
{
    BaseExp::init(input, gwstate, nullptr, nullptr, logging, _proof.depth() + (security() ? 1 : 0));
    _state_update_period = 1;
    _logging->set_prefix(input->display_text() + " ");
    if (security())
    {
        _rnd_seed = _security_seed;
        _rnd_seed.arithmetic().alloc(_rnd_seed, _rnd_seed.size() + 2);
        *(double *)(_rnd_seed.data() + _rnd_seed.size()) = getHighResTimer();
        _rnd_seed.arithmetic().init(_rnd_seed.data(), _rnd_seed.size() + 2, _rnd_seed);
        _logging->info("random seed: %s.\n", _rnd_seed.to_string().data());
    }
}

void ProofBuild::execute()
{
    int M, t, i, j, k;
    Proof::Product state_d;
    GWNum X(gw());
    GWNum Y(gw());
    GWNum D(gw());
    GWNum T(gw());
    Giant exp;
    int len = 0;
    Giant a_power;
    std::vector<Giant> tree;
    std::vector<Giant> h;

    t = _proof.depth();
    X = _proof.r_0();
    if (_proof.Li())
    {
        len = _proof.r_exp().bitlen() - 1;
        _proof.r_exp().arithmetic().substr(_proof.r_exp(), len - _proof.points()[1], _proof.points()[1], a_power);
    }

    if (state() == nullptr)
        _state.reset(new State(0, std::move(_proof.r_count())));
    Y = state()->X();

    M = _proof.points()[_proof.count()];
    for (i = 0; i < t; i++, commit_execute<State>(i, Y), M >>= 1)
    {
        _proof.read_product(i, state_d, *_logging);
        D = state_d.X();
        h.emplace_back(GiantsArithmetic::default_arithmetic(), 4);
        hash_giants(_gwstate->fingerprint, state()->X(), state_d.X(), h[i]);
        make_prime(h[i]);

        if (_proof.Li())
        {
            a_power *= h[i];
            while (tree.size() < i)
                tree.emplace_back();
            for (j = 0; j < (1 << i); j++)
            {
                k = (1 + j*2) << (t - i - 1);
                _proof.r_exp().arithmetic().substr(_proof.r_exp(), len - _proof.points()[k + 1], _proof.points()[k + 1] - _proof.points()[k], exp);
                for (k = 1; k <= i; k++)
                    if ((j & (1 << (k - 1))) == 0)
                    {
                        tree[i - k] = exp;
                        break;
                    }
                    else
                        exp += tree[i - k]*h[i - k];
            }
            a_power += exp;
        }

        exp = h[i];
        if (M%2 != 0 && !_proof.Li())
            exp *= _input->gb();
        if (M%2 != 0 && _proof.Li())
            exp *= 2;
        exp_gw(gw().carefully(), exp, X, T = X, 0);
        gw().carefully().mul(D, X, X, 0);

        exp_gw(gw().carefully(), h[i], D, T = D, 0);
        gw().carefully().mul(D, Y, Y, 0);
    }

    D = _proof.r_0();
    if (!_rnd_seed.empty())
    {
        _raw_res64 = state()->X().to_res64();

        exp.arithmetic().rnd_seed(_rnd_seed);
        exp.arithmetic().rnd(exp, 64);
        make_prime(exp, 1000000);

        exp_gw(gw().carefully(), exp, X, T = X, 0);
        exp_gw(gw().carefully(), exp, Y, T = Y, 0);
        if (_proof.Li())
            exp_gw(gw().carefully(), exp, D, T = D, 0);

        commit_execute<State>(t + 1, Y);
    }

    if (state()->X() == 0)
    {
        _logging->error("invalid proof, the certificate is zero.\n");
        throw TaskAbortException();
    }

    Proof::Certificate cert;
    if (!_proof.Li())
        cert.set(M, X);
    else
        cert.set(M, X, a_power, D);
    _proof.file_cert()->write(cert);
    _proof.file_cert()->free_buffer();

    done();
}
