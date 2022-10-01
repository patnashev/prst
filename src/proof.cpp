
#include <cmath>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "proof.h"
#include "md5.h"
#include "integer.h"

using namespace arithmetic;

Proof::Proof(int op, int count, InputNum& input, Params& params, Logging& logging) : _op(op), _count(count)
{
    if ((op == SAVE || op == BUILD) && (count & (count - 1)) != 0)
    {
        logging.error("proof count is not a power of 2.\n");
        exit(0);
    }

    if (op == SAVE)
        _task.reset(new ProofSave(*this));
    if (op == BUILD)
        _task.reset(new ProofBuild(*this, params.ProofSecuritySeed));
    if (op == CERT)
    {
        _M = count;
        _points.push_back(count);
        bool CheckGerbicz = params.CheckGerbicz ? params.CheckGerbicz.value() : input.b() == 2;

        if (!CheckGerbicz)
        {
            MultipointExp* task;
            _task.reset(task = new MultipointExp(input.gb(), _points, nullptr));
            if (params.SlidingWindow)
                task->_W = params.SlidingWindow.value();
            logging.progress().add_stage(task->cost());
        }
        else
        {
            int checks = params.GerbiczCount ? params.GerbiczCount.value() : 1;
            GerbiczCheckExp* task;
            _task.reset(task = params.GerbiczL ? new GerbiczCheckExp(input.gb(), _points[0], checks, params.GerbiczL.value()) : new GerbiczCheckExp(input.gb(), _points[0], checks));
            if (params.SlidingWindow)
                task->_W = params.SlidingWindow.value();
            logging.progress().add_stage(task->cost());
        }
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
        _taskRoot.reset(new SlowExp(std::move(exp)));
    }

    if (_task)
        _task->set_error_check(params.CheckNear && params.CheckNear.value(), !params.Check || params.Check.value());
    if (_taskRoot)
        _taskRoot->set_error_check(params.CheckNear && params.CheckNear.value(), !params.Check || params.Check.value());
}

int Proof::read_cert_power(File& file_cert)
{
    std::unique_ptr<Reader> reader(file_cert.get_reader());
    if (!reader)
        return 0;
    if (reader->type() != Certificate::TYPE)
        return 0;
    int power = 0;
    reader->read(power);
    return power;
}

void Proof::calc_points(int iterations, InputNum& input, Params& params, Logging& logging)
{
    if (params.GerbiczCount)
        if (params.GerbiczCount.value() > _count)
            params.ProofChecksPerPoint = params.GerbiczCount.value()/_count;
        else
        {
            _points_per_check = _count/params.GerbiczCount.value();
            params.ProofPointsPerCheck = _points_per_check;
        }
    if ((input.b() != 2 || input.c() != 1) && _points_per_check > 1)
    {
        int i;
        int iters = iterations*_points_per_check/_count;
        if (!params.GerbiczL)
        {
            int L = _points_per_check*(int)sqrt(iters);
            int L2 = iters - iters%L;
            for (i = L + _points_per_check; i*i < 2*iters*_points_per_check*_points_per_check; i += _points_per_check)
                if (L2 < iters - iters%i)
                {
                    L = i;
                    L2 = iters - iters%i;
                }
            params.GerbiczL = L/_points_per_check;
            _M = L2/_points_per_check;
        }
        else
            _M = (iters - iters%(params.GerbiczL.value()*_points_per_check))/_points_per_check;
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

        read_point(_count, *state, logging);
        _r_count = state->X();

        task->set_error_check(false, true);
        task->init_state(state.release());
        return;
    }

    logging.info("Saving %d proof points.\n", count());
    GerbiczCheckMultipointExp* taskGerbiczCheck = dynamic_cast<GerbiczCheckMultipointExp*>(task);
    int point = _count;
    while (point >= 0)
    {
        if (taskGerbiczCheck != nullptr)
            point -= point%_points_per_check;
        if (_file_points[point]->read(*state) && state->iteration() == _points[point])
        {
            if (task->state() == nullptr || task->state()->iteration() < state->iteration())
                task->init_state(state.release());
            return;
        }
        point--;
    }
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

void Proof::on_point(int index, arithmetic::Giant& X)
{
    if (index > _count)
        return;
    BaseExp::State state(_points[index], std::move(X));
    _file_points[index]->write(state);
    if (!_cache_points)
        _file_points[index]->free_buffer();
    X = std::move(state.X());
}

void Proof::run(InputNum& input, arithmetic::GWState& gwstate, File& file_cert, File& file_checkpoint, File& file_recoverypoint, Logging& logging)
{
    logging.info("Verifying certificate of %s, %d iterations.\n", input.display_text().data(), _points[0]);

    MultipointExp* taskMultipoint = dynamic_cast<MultipointExp*>(_task.get());
    if (taskMultipoint != nullptr)
    {
        GerbiczCheckMultipointExp* taskGerbiczCheck = dynamic_cast<GerbiczCheckMultipointExp*>(taskMultipoint);
        if (taskGerbiczCheck != nullptr)
            taskGerbiczCheck->init(&input, &gwstate, &file_checkpoint, &file_recoverypoint, &logging);
        else
            taskMultipoint->init(&input, &gwstate, &file_checkpoint, &logging);
        if (taskMultipoint->state() == nullptr)
        {
            std::unique_ptr<Certificate> cert(read_state<Certificate>(&file_cert));
            file_cert.free_buffer();
            taskMultipoint->init_state(new BaseExp::State(0, std::move(cert->X())));
        }
    }
    _task->run();

    _res64 = _task->state()->X().to_res64();
    logging.result(false, "%s certificate RES64: : %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), _task->timer());
    logging.result_save(input.input_text() + " certificate RES64: " + _res64 + ", time: " + std::to_string((int)_task->timer()) + " s.\n");

    logging.progress().next_stage();
    file_checkpoint.clear();
    file_recoverypoint.clear();
}

void Proof::run(InputNum& input, arithmetic::GWState& gwstate, Logging& logging, Giant* X)
{
    SlowExp* taskRoot = dynamic_cast<SlowExp*>(_taskRoot.get());
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
    int i, M, t;
    Proof::Product state_d;
    GWNum X(gw());
    GWNum Y(gw());
    GWNum D(gw());
    GWNum T(gw());
    Giant h(GiantsArithmetic::default_arithmetic(), 4);
    Giant exp;

    t = _proof.depth();
    X = _proof.r_0();

    if (state() == nullptr)
        _state.reset(new State(0, std::move(_proof.r_count())));
    Y = state()->X();

    M = _proof.points()[_proof.count()];
    for (i = 0; i < t; i++, commit_execute<State>(i, Y), M >>= 1)
    {
        _proof.read_product(i, state_d, *_logging);
        D = state_d.X();
        hash_giants(_gwstate->fingerprint, state()->X(), state_d.X(), h);
        make_prime(h);

        exp = h;
        if (M%2 != 0)
            exp *= _input->gb();
        exp_gw(gw().carefully(), exp, X, T = X, 0);
        gw().carefully().mul(D, X, X, 0);

        exp_gw(gw().carefully(), h, D, T = D, 0);
        gw().carefully().mul(D, Y, Y, 0);
    }

    if (!_rnd_seed.empty())
    {
        _raw_res64 = state()->X().to_res64();

        exp.arithmetic().rnd_seed(_rnd_seed);
        exp.arithmetic().rnd(exp, 64);
        make_prime(exp, 1000000);

        exp_gw(gw().carefully(), exp, X, T = X, 0);
        exp_gw(gw().carefully(), exp, Y, T = Y, 0);

        commit_execute<State>(t + 1, Y);
    }

    if (state()->X() == 0)
    {
        _logging->error("invalid proof, the certificate is zero.\n");
        throw TaskAbortException();
    }

    Proof::Certificate cert;
    cert.set(M, X);
    _proof.file_cert()->write(cert);
    _proof.file_cert()->free_buffer();

    done();
}
