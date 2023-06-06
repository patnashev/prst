
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
        _task.reset(new ProofSave());
    if (op == BUILD)
        _task.reset(new ProofBuild(params.ProofSecuritySeed));
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
                _taskA.reset(new SlidingWindowExp(std::move(cert.a_power())));
            else
                _taskA.reset(new LiCheckExp(std::move(cert.a_power()), checks, params.StrongL ? params.StrongL.value() : 0));
            if (params.SlidingWindow)
                _taskA->_W = params.SlidingWindow.value();
            logging.progress().add_stage(_taskA->cost());
        }
        logging.progress().add_stage(task->cost());
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
    if (params.CheckStrong && params.CheckStrong.value() && params.StrongCount)
        if (params.StrongCount.value() > _count)
            params.ProofChecksPerPoint = params.StrongCount.value()/_count;
        else
            params.ProofPointsPerCheck = _count/params.StrongCount.value();
    int points_per_check = params.ProofPointsPerCheck ? params.ProofPointsPerCheck.value() : 1;
    bool value_points = input.type() != input.KBNC || input.k() == 0 || input.b() == 0 || input.b() == 2;

    /*if (input.b() != 2 && params.ProofPointsPerCheck && !Li())
    {
        int i;
        int iters = iterations*points_per_check/_count;
        if (!params.StrongL)
        {
            int L = points_per_check*(int)std::sqrt(iters);
            int L2 = iters - iters%L;
            for (i = L + points_per_check; i*i < 2*iters*points_per_check*points_per_check; i += points_per_check)
                if (L2 < iters - iters%i)
                {
                    L = i;
                    L2 = iters - iters%i;
                }
            params.StrongL = L/points_per_check;
            _M = L2/points_per_check;
        }
        else
            _M = (iters - iters%(params.StrongL.value()*points_per_check))/points_per_check;
        _points.reserve(_count + 2);
        for (i = 0; i <= _count; i++)
            _points.emplace_back(i*_M, i%points_per_check == 0 || i == _count, value_points || i == _count);
    }
    else*/
    {
        _points.reserve(_count + 1);
        _points.emplace_back(0);
        for (int i = 1; i < _count; i++)
        {
            _M = iterations;
            _points.emplace_back(0, i%points_per_check == 0, value_points);
            for (int j = _count/2; j > 0 && (i & (j*2 - 1)) != 0; j >>= 1)
            {
                _M /= 2;
                if ((i & j) != 0)
                    _points.back().pos += _M;
                if ((iterations & (_count/j/2)) != 0)
                    _points.back().pos++;
            }
        }
        _points.push_back(iterations);
    }
    logging.report_param("M", _M);
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
    std::unique_ptr<BaseExp::StateValue> state(new BaseExp::StateValue());
    std::unique_ptr<BaseExp::StateSerialized> state_serialized(new BaseExp::StateSerialized());

    if (op() == BUILD)
    {
        logging.info("Building certificate from %d products.\n", depth());
        logging.set_prefix(input.display_text() + " ");

        if (!Li())
        {
            read_point(0, *state, logging);
            double timer = getHighResTimer();
            Giant tmp;
            tmp = a;
            tmp.arithmetic().powermod(tmp, input.gk(), *gwstate.N, tmp);
            if (tmp != state->value())
            {
                logging.error("invalid a^k.\n");
                throw TaskAbortException();
            }
            timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
            logging.info("a^k is correct. Time: %.1f s.\n", timer);
            _r_0 = std::move(state->value());
        }
        else
        {
            _r_0 = a;
            _r_exp = &task->exp();
        }

        read_point(_count, *state, logging);
        _r_count = state->value();

        task->set_error_check(false, true);
        task->init_state(state.release());
        return;
    }

    logging.info("Saving %d proof points.\n", count());
    int point = _count;
    while (point >= 0)
    {
        if (point == 0 && Li())
        {
            return;
        }
        BaseExp::State* tmp_state;
        if (_points[point].check && (tmp_state = BaseExp::State::read_file(_file_points[point], state.get(), state_serialized.get())) != nullptr && tmp_state->iteration() == _points[point].pos)
        {
            if (task->state() == nullptr || task->state()->iteration() < tmp_state->iteration())
            {
                task->init_state(tmp_state);
                if (tmp_state == state.get())
                    state.release();
                if (tmp_state == state_serialized.get())
                    state_serialized.release();
            }
            return;
        }
        point--;
    }
    task->init_state(nullptr);
}

void Proof::read_point(int index, TaskState& state, Logging& logging)
{
    if (!_file_points[index]->read(state) || state.iteration() != _points[index].pos)
    {
        logging.error("%s is missing or corrupt.\n", _file_points[index]->filename().data());
        throw TaskAbortException();
    }
    _file_points[index]->free_buffer();
}

BaseExp::State* Proof::read_point(int index, BaseExp::StateValue* state_value, BaseExp::StateSerialized* state_serialized, Logging& logging)
{
    BaseExp::State* state = BaseExp::State::read_file(_file_points[index], state_value, state_serialized);
    if (state == nullptr || state->iteration() != _points[index].pos)
    {
        logging.error("%s is missing or corrupt.\n", _file_points[index]->filename().data());
        throw TaskAbortException();
    }
    _file_points[index]->free_buffer();
    return state;
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

bool Proof::on_point(int index, BaseExp::State* state)
{
    if (index > _count)
        return false;
    GWASSERT(index < _count || dynamic_cast<BaseExp::StateValue*>(state) != nullptr);
    _file_points[index]->write(*state);
    if (!_cache_points)
        _file_points[index]->free_buffer();
    return true;
}

void Proof::run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging)
{
    MultipointExp* task = dynamic_cast<MultipointExp*>(_task.get());
    if (task == nullptr)
        return;

    double timer = 0;
    Giant tail;
    if (Li())
    {
        File* file_checkpoint_a(file_checkpoint.add_child("a", file_checkpoint.fingerprint()));
        File* file_recoverypoint_a(file_recoverypoint.add_child("a", file_recoverypoint.fingerprint()));
        logging.info("Verifying certificate of %s, complexity = %d+%d.\n", input.display_text().data(), (int)logging.progress().costs()[0], (int)logging.progress().costs()[1]);
        logging.set_prefix(input.display_text() + " ");

        LiCheckExp* taskACheck = dynamic_cast<LiCheckExp*>(_taskA.get());
        if (taskACheck != nullptr)
            taskACheck->init(&input, &gwstate, file_checkpoint_a, file_recoverypoint_a, &logging, std::move(_r_0));
        else
            _taskA->init(&input, &gwstate, file_checkpoint_a, &logging, std::move(_r_0));
        _taskA->run();
        timer = _taskA->timer();
        tail = std::move(*_taskA->result());

        logging.progress().next_stage();
        file_checkpoint_a->clear();
        file_recoverypoint_a->clear();
    }
    else
    {
        logging.info("Verifying certificate of %s, complexity = %d.\n", input.display_text().data(), (int)logging.progress().cost_total());
        logging.set_prefix(input.display_text() + " ");
    }

    GerbiczCheckExp* taskCheck = dynamic_cast<GerbiczCheckExp*>(_task.get());
    if (taskCheck != nullptr)
        taskCheck->init(&input, &gwstate, &file_checkpoint, &file_recoverypoint, &logging, std::move(tail));
    else
        task->init_smooth(&input, &gwstate, &file_checkpoint, &logging, std::move(tail));
    if (task->state() == nullptr)
    {
        task->init_state(new BaseExp::StateValue(0, std::move(_r_count)));
        task->state()->set_written();
    }
    _task->run();
    timer += _task->timer();
    logging.set_prefix("");

    _res64 = task->result()->to_res64();
    logging.result(false, "%s certificate RES64: %s, time: %.1f s.\n", input.display_text().data(), _res64.data(), timer);
    logging.result_save(input.input_text() + " certificate RES64: " + _res64 + ", time: " + std::to_string((int)timer) + " s.\n");

    logging.progress().next_stage();
    file_checkpoint.clear();
    file_recoverypoint.clear();
}

void Proof::run(InputNum& input, arithmetic::GWState& gwstate, Logging& logging, Giant* X)
{
    if (_taskRoot && X != nullptr)
    {
        _taskRoot->init(&input, &gwstate, &logging, std::move(*X));
        _taskRoot->run();
        if (*_taskRoot->result() == 1)
        {
            logging.error("%s roots of unity check failed.\n", input.display_text().data());
            throw TaskAbortException();
        }
    }
    ProofSave* taskSave = dynamic_cast<ProofSave*>(_task.get());
    if (taskSave != nullptr)
    {
        logging.set_prefix(input.display_text() + " ");
        taskSave->init(&input, &gwstate, &logging, this);
        taskSave->run();
        logging.set_prefix("");
        logging.progress().next_stage();

        _res64 = taskSave->state()->Y().to_res64();
        logging.info("%s compressed %d points to %d products, time: %.1f s.\n", input.display_text().data(), _count, _file_products.size(), _task->timer());
        logging.result(false, "%s raw certificate RES64: %s.\n", input.display_text().data(), _res64.data());
        logging.result_save(input.input_text() + " raw certificate RES64: " + _res64 + ", time: " + std::to_string((int)_task->timer()) + " s.\n");
    }
    ProofBuild* taskBuild = dynamic_cast<ProofBuild*>(_task.get());
    if (taskBuild != nullptr)
    {
        logging.set_prefix(input.display_text() + " ");
        taskBuild->init(&input, &gwstate, &logging, this);
        taskBuild->run();
        logging.set_prefix("");
        logging.progress().next_stage();

        if (taskBuild->security())
        {
            _res64 = taskBuild->raw_res64();
            logging.result(false, "%s raw certificate RES64: %s.\n", input.display_text().data(), _res64.data());
            logging.result_save(input.input_text() + " raw certificate RES64: " + _res64 + ".\n");
        }
        _res64 = taskBuild->state()->Y().to_res64();
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

void ProofSave::init(InputNum* input, arithmetic::GWState* gwstate, Logging* logging, Proof* proof)
{
    InputTask::init(input, gwstate, nullptr, nullptr, logging, proof->count());
    _state_update_period = 0;
    _logging->set_prefix(input->display_text() + " ");
    _proof = proof;
}

void ProofSave::done()
{
    InputTask::done();
    _logging->set_prefix("");
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
    Proof::Product product;
    BaseExp::StateValue state_v;
    BaseExp::StateSerialized state_s;
    GWNum Y(gw());
    GWNum D(gw());
    GWNum T(gw());
    std::vector<GWNum> tree;
    std::vector<Giant> h;

    t = _proof->depth();
    tree.reserve(t);
    h.reserve(t);

    if (state() == nullptr)
    {
        _proof->read_point(_proof->count(), state_v, *_logging);
        _state.reset(new Proof::State(0, std::move(state_v.value())));
    }
    Y = state()->Y();
    h = state()->h();

    for (i = state()->iteration(); i < t; i++, commit_execute<Proof::State>(i, Y, h))
    {
        if (_proof->file_products()[i]->read(product))
        {
            _proof->file_products()[i]->free_buffer();
            D = product.value();
        }
        else
        {
            if (i == 0)
            {
                _proof->read_point(_proof->count()/2, &state_v, &state_s, *_logging)->to_GWNum(D);
            }
            else
            {
                while (tree.size() < i)
                    tree.emplace_back(gw());
                for (j = 0; j < (1 << i); j++)
                {
                    k = (1 + j*2) << (t - i - 1);
                    _proof->read_point(k, &state_v, &state_s, *_logging)->to_GWNum(D);

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
                check();
            }
            product.set(i, D);

            _proof->file_products()[i]->write(product);
            _proof->file_products()[i]->free_buffer();
            on_state();
        }

        h.emplace_back(GiantsArithmetic::default_arithmetic(), 4);
        hash_giants(_gwstate->fingerprint, state()->Y(), product.value(), h[i]);
        make_prime(h[i]);

        gw().fft(D, D);
        exp_gw(gw(), h[i], D, T = D, GWMUL_STARTNEXTFFT);
        gw().mul(D, Y, Y, 0);
    }
    tree.clear();

    if (_input->need_mod())
        _input->mod(state()->Y(), state()->Y());
    done();
}

void ProofBuild::init(InputNum* input, arithmetic::GWState* gwstate, Logging* logging, Proof* proof)
{
    InputTask::init(input, gwstate, nullptr, nullptr, logging, proof->depth() + (security() ? 1 : 0));
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
    _proof = proof;
}

void ProofBuild::done()
{
    InputTask::done();
    _logging->set_prefix("");
}

void ProofBuild::execute()
{
    int M, t, i, j, k;
    Proof::Product product;
    GWNum X(gw());
    GWNum Y(gw());
    GWNum D(gw());
    GWNum T(gw());
    Giant exp;
    int len = 0;
    Giant a_power;
    std::vector<Giant> tree;
    std::vector<Giant> h;

    t = _proof->depth();
    if (_proof->Li())
    {
        len = _proof->r_exp().bitlen() - 1;
        _proof->r_exp().arithmetic().substr(_proof->r_exp(), len - _proof->points()[1].pos, _proof->points()[1].pos, a_power);
    }

    if (state() == nullptr)
    {
        X = _proof->r_0();
        _state.reset(new Proof::State(0, X, std::move(_proof->r_count()), std::move(a_power)));
    }
    else
        X = state()->X();
    Y = state()->Y();
    a_power = state()->exp();
    h = state()->h();

    M = _proof->points()[_proof->count()].pos >> state()->iteration();
    for (i = state()->iteration(); i < t; i++, commit_execute<Proof::State>(i, X, Y, a_power, h), M >>= 1)
    {
        _proof->read_product(i, product, *_logging);
        D = product.value();
        h.emplace_back(GiantsArithmetic::default_arithmetic(), 4);
        hash_giants(_gwstate->fingerprint, state()->Y(), product.value(), h[i]);
        make_prime(h[i]);

        if (_proof->Li())
        {
            a_power *= h[i];
            while (tree.size() < i)
                tree.emplace_back();
            for (j = 0; j < (1 << i); j++)
            {
                k = (1 + j*2) << (t - i - 1);
                _proof->r_exp().arithmetic().substr(_proof->r_exp(), len - _proof->points()[k + 1].pos, _proof->points()[k + 1].pos - _proof->points()[k].pos, exp);
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
        if (M%2 != 0 && !_proof->Li())
            exp *= _input->gb();
        if (M%2 != 0 && _proof->Li())
            exp *= 2;
        exp_gw(gw().carefully(), exp, X, T = X, 0);
        gw().carefully().mul(D, X, X, 0);

        exp_gw(gw().carefully(), h[i], D, T = D, 0);
        gw().carefully().mul(D, Y, Y, 0);
    }
    if (_input->need_mod())
        _input->mod(state()->Y(), state()->Y());

    D = _proof->r_0();
    if (!_rnd_seed.empty())
    {
        _raw_res64 = state()->Y().to_res64();

        exp.arithmetic().rnd_seed(_rnd_seed);
        exp.arithmetic().rnd(exp, 64);
        make_prime(exp, 1000000);

        exp_gw(gw().carefully(), exp, X, T = X, 0);
        exp_gw(gw().carefully(), exp, Y, T = Y, 0);
        if (_proof->Li())
            exp_gw(gw().carefully(), exp, D, T = D, 0);

        commit_execute<Proof::State>(t + 1, X, Y, a_power, h);
    }
    if (_input->need_mod())
        _input->mod(state()->Y(), state()->Y());

    if (state()->Y() == 0)
    {
        _logging->error("invalid proof, the certificate is zero.\n");
        throw TaskAbortException();
    }

    Proof::Certificate cert;
    if (!_proof->Li())
        cert.set(M, X);
    else
        cert.set(M, X, a_power, D);
    _proof->file_cert()->write(cert);
    _proof->file_cert()->free_buffer();

    done();
}
