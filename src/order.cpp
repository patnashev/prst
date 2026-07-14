
#include <algorithm>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "integer.h"

#include "order.h"

using namespace arithmetic;

Order::Order(InputNum& input, Options& options, Logging& logging) : Run("Order", input, options)
{
    _factors = input.factors();
    Giant ga = options.OrderA->value();
    create_tasks(ga, logging, false);

    if (_task && ga <= GWMULBYCONST_MAX)
        logging.report_param("maxmulbyconst", (int)ga.data()[0]);

    _fingerprint = File::unique_fingerprint(input.fingerprint(), std::to_string(options.OrderA->fingerprint()));
}

MultipointExp* Order::create_smooth_task(arithmetic::Giant& base, int power)
{
    bool CheckStrong = _options.CheckStrong ? _options.CheckStrong.value() : true;
    int checks = _options.StrongCount ? _options.StrongCount.value() : 16;

    MultipointExp* task;
    if (CheckStrong)
    {
        task = new GerbiczCheckExp(base, power, checks, std::bind(&Order::on_point, this, std::placeholders::_1, std::placeholders::_2));
        for (auto& p : task->points())
            p.value = true;
    }
    else
    {
        int len = power/checks;
        if (len == 0)
            len = 1;
        std::vector<MultipointExp::Point> points;
        for (int i = 0; i <= checks && len*i <= power; i++)
            points.emplace_back(len*i);
        if (points.back().pos != power)
            points.emplace_back(power);
        task = new MultipointExp(base, true, points, std::bind(&Order::on_point, this, std::placeholders::_1, std::placeholders::_2));
    }
    return task;
}

void Order::create_tasks(Giant& a, Logging& logging, bool restart)
{
    bool CheckStrong = _options.CheckStrong ? _options.CheckStrong.value() : true;
    int checks = _options.StrongCount ? _options.StrongCount.value() : 16;

    _tasks_smooth.clear();
    Giant exp;
    exp = 1;
    _task_exp_str.clear();
    for (auto& factor : _factors)
    {
        if (factor.second <= _sub)
            continue;
        if (factor.first == 2 || restart)
        {
            _tasks_smooth.emplace_back(create_smooth_task(factor.first, factor.second - _sub));
        }
        else
        {
            exp *= power(factor.first, factor.second - _sub);
            if (!_task_exp_str.empty())
                _task_exp_str += "*";
            _task_exp_str += factor.first.to_string();
            _task_exp_str += "^";
            _task_exp_str += std::to_string(factor.second - _sub);
        }
    }
    for (auto& factor : _order)
    {
        exp *= power(factor.first, factor.second);
        if (!_task_exp_str.empty())
            _task_exp_str += "*";
        _task_exp_str += factor.first.to_string();
        _task_exp_str += "^";
        _task_exp_str += std::to_string(factor.second);
    }
    if (exp == 1)
        _task.reset();
    else
    {
        if (CheckStrong && a <= GWMULBYCONST_MAX)
            _task.reset(new FastLiCheckExp(exp, exp.bitlen() > 100 ? checks : 1));
        else if (CheckStrong)
            _task.reset(new LiCheckExp(exp, exp.bitlen() > 100 ? checks : 1));
        else if (a <= GWMULBYCONST_MAX)
            _task.reset(new FastExp(exp));
        else
            _task.reset(new SlidingWindowExp(exp));
        logging.progress().add_stage(_task->cost());
    }
    for (auto& task_smooth : _tasks_smooth)
        logging.progress().add_stage(task_smooth->cost());

    _tasks.clear();
    exp = 1;
    for (auto& factor : _factors)
    {
        int n = (factor.second <= _sub ? factor.second : _sub);
        exp *= power(factor.first, n);
        _tasks.emplace_back(factor.first, factor.second - n, n);
    }
    _task_check.reset(new CarefulExp(exp));
    for (auto& factor : _tasks)
    {
        factor.task_sub.reset(new CarefulExp(exp/power(factor.b, factor.n)));
        factor.task_factor.reset(new CarefulExp(factor.b));
    }
}

bool Order::on_point(int index, BaseExp::State* state)
{
    if (dynamic_cast<BaseExp::StateValue*>(state)->value() == 1)
    {
        _task_break = index;
        throw TaskAbortException();
    }
    return false;
}

void Order::run(arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging)
{
    Giant ga = _options.OrderA->value();
    logging.info("Computing multiplicative order of %s modulo prime %s.\n", _options.OrderA->display_text().data(), input.display_text().data());
    if (gwstate.information_only)
        throw TaskAbortException();
    logging.set_prefix("ord(" + _options.OrderA->display_text() + ") mod " + input.display_text() + " ");

    while (!_factors.empty())
    {
        Giant sub_val;
        sub_val = ga;

        if (_task)
        {
            logging.info("raising to power %s.\n", _task_exp_str.data());
            if (FastExp* task = dynamic_cast<FastExp*>(_task.get()))
                task->init(&input, &gwstate, nullptr, &logging, ga.data()[0]);
            if (SlidingWindowExp* task = dynamic_cast<SlidingWindowExp*>(_task.get()))
                task->init(&input, &gwstate, nullptr, &logging, ga);
            if (FastLiCheckExp* task = dynamic_cast<FastLiCheckExp*>(_task.get()))
                task->init(&input, &gwstate, nullptr, nullptr, &logging, ga.data()[0]);
            else if (LiCheckExp* task = dynamic_cast<LiCheckExp*>(_task.get()))
                    task->init(&input, &gwstate, nullptr, nullptr, &logging, ga);

            _task->run();
            sub_val = *_task->result();
            logging.progress().next_stage();
        }
        if (sub_val == 1 && !_order.empty())
        {
            _factors.clear();
            break;
        }
        for (auto& task_smooth : _tasks_smooth)
        {
            if (sub_val == 1)
            {
                auto it = _factors.begin();
                for (; it != _factors.end() && it->first != task_smooth->exp(); it++);
                _factors.erase(it);
                logging.progress().next_stage();
                continue;
            }
            logging.info("raising to power %s^%d.\n", task_smooth->b().to_string().data(), task_smooth->points().back());
            if (GerbiczCheckExp* task = dynamic_cast<GerbiczCheckExp*>(task_smooth.get()))
                task->init(&input, &gwstate, nullptr, nullptr, &logging);
            else
                task_smooth->init_smooth(&input, &gwstate, nullptr, &logging);
            if (task_smooth->state() == nullptr)
                task_smooth->init_state(new BaseExp::StateValue(0, sub_val));

            try
            {
                _task_break = -1;
                task_smooth->run();
                sub_val = *task_smooth->result();
            }
            catch (const TaskAbortException&)
            {
                if (_task_break == -1)
                    throw;
                sub_val = 1;
                auto it = _factors.begin();
                for (; it != _factors.end() && it->first != task_smooth->exp(); it++);
                it->second = task_smooth->points()[_task_break].pos + _sub;
            }
            logging.progress().next_stage();
        }

        if (sub_val == 1)
        {
            for (auto it = _factors.begin(); it != _factors.end(); it++)
                if (it->second > _sub)
                    it->second -= _sub;
            create_tasks(ga, logging, true);
            continue;
        }

        logging.info("computing order for each factor.\n");

        _task_check->init_giant(&input, &gwstate, &logging, std::move(sub_val));
        _task_check->run();
        sub_val = std::move(_task_check->X0());
        if (*_task_check->result() != 1)
        {
            logging.set_prefix("");
            logging.error("%s is not prime.\n", input.display_text().data());
            throw TaskAbortException();
        }

        for (auto& factor : _tasks)
        {
            auto it = _factors.begin();
            for (; it != _factors.end() && it->first != factor.b; it++);

            factor.task_sub->init_giant(&input, &gwstate, &logging, std::move(sub_val));
            factor.task_sub->run();
            sub_val = std::move(factor.task_sub->X0());

            Giant cur_val = std::move(*factor.task_sub->result());
            if (cur_val == 1)
            {
                if (it->second <= _sub)
                    _factors.erase(it);
                else
                    it->second -= _sub;
                continue;
            }

            int ord = factor.ord;
            for (int i = 0; i < factor.n && cur_val != 1; i++, ord++)
            {
                factor.task_factor->init_giant(&input, &gwstate, &logging, std::move(cur_val));
                factor.task_factor->run();
                cur_val = std::move(*factor.task_factor->result());
            }

            if (cur_val != 1)
            {
                logging.set_prefix("");
                logging.error("%s is not prime.\n", input.display_text().data());
                throw TaskAbortException();
            }

            _order.emplace_back(factor.b, ord);
            _factors.erase(it);
        }

        if (!_factors.empty())
            create_tasks(ga, logging, true);
    }

    Giant order_div;
    order_div = 1;
    for (auto& factor : _order)
    {
        auto it = input.factors().begin();
        for (; it != input.factors().end() && it->first != factor.first; it++);
        if (it->second != factor.second)
            order_div *= power(factor.first, it->second - factor.second);

        if (!_result.empty())
            _result += "*";
        _result += factor.first.to_string();
        if (factor.second > 1)
        {
            _result += "^";
            _result += std::to_string(factor.second);
        }
    }
    for (auto& factor : input.factors())
    {
        auto it = _order.begin();
        for (; it != _order.end() && it->first != factor.first; it++);
        if (it == _order.end())
            order_div *= power(factor.first, factor.second);
    }
    if (order_div.size() == 1)
    {
        _result += " = (N-1)/";
        _result += order_div.to_string();
    }

    logging.set_prefix("");
    logging.result(true, "ord(%s) mod %s = %s.\n", _options.OrderA->display_text().data(), input.display_text().data(), _result.data());
    logging.result_save("ord(" + _options.OrderA->display_text() + ") mod " + input.input_text() + " = " + _result + ".\n");

    file_checkpoint.clear(true);
    file_recoverypoint.clear(true);
}

FermatDivisor::FermatDivisor(InputNum& input, Options& options, Logging& logging) : Order("FermatDivisor", input, options)
{
    GWASSERT(input.factors()[0].first == 2);
    
    int limit = options.Divides == "F" ? 2 : options.DividesLimit ? options.DividesLimit.value() : 12;
    for (int i = 2; i <= limit; i++)
    {
        _bases.emplace_back(i);
        if (is_prime(i))
            create_task(logging, _bases.back());
    }

    _fingerprint = File::unique_fingerprint(input.fingerprint(), "div");
}

void FermatDivisor::create_task(Logging& logging, Base& base)
{
    base.task.reset();
    base.exp = 1;

    int n = input.factors()[0].second;
    if (logging.progress().param_int("base_" + base.str) != 0)
        n = logging.progress().param_int("base_" + base.str);
    if (n <= _sub)
    {
        base.power = 0;
        base.exp <<= n;
        base.sub_val = base.base;
    }
    else
    {
        base.power = n - _sub;
        base.exp <<= _sub;
        base.task.reset(create_smooth_task(input.factors()[0].first, base.power));
        logging.progress().add_stage(base.task->cost());
    }
}

void FermatDivisor::run_task(arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Base& base)
{
    while (base.task)
    {
        File* checkpoint = file_checkpoint.add_child(base.str, File::unique_fingerprint(file_checkpoint.fingerprint(), base.str + "." + std::to_string(base.power)));
        File* recoverypoint = file_recoverypoint.add_child(base.str, File::unique_fingerprint(file_recoverypoint.fingerprint(), base.str + "." + std::to_string(base.power)));

        logging.info("raising %s to power 2^%d.\n", base.str.data(), base.power);
        if (GerbiczCheckExp* taskCheck = dynamic_cast<GerbiczCheckExp*>(base.task.get()))
            taskCheck->init(&input, &gwstate, checkpoint, recoverypoint, &logging);
        else
            base.task->init_smooth(&input, &gwstate, checkpoint, &logging);
        if (base.task->state() == nullptr)
            base.task->init_state(new BaseExp::StateValue(0, base.base));

        try
        {
            _task_break = -1;
            base.task->run();
            base.task->write_state();
            base.sub_val = std::move(*base.task->result());

            base.task.reset();
            checkpoint->free_buffer();
            recoverypoint->free_buffer();
        }
        catch (const TaskAbortException&)
        {
            if (_task_break == -1)
                throw;
            logging.report_param("base_" + base.str, base.task->points()[_task_break].pos);
            logging.progress_save();
            checkpoint->clear();
            recoverypoint->clear();
            
            create_task(logging, base);
            logging.progress().costs()[logging.progress().cur_stage()] *= logging.progress().progress_stage();
            if (base.task)
                std::swap(logging.progress().costs()[logging.progress().cur_stage() + 1], logging.progress().costs().back());
        }
        logging.progress().next_stage();
    }
}

void FermatDivisor::run(arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging)
{
    int limit = _options.Divides == "F" ? 2 : _options.DividesLimit ? _options.DividesLimit.value() : 12;
    if (limit > 2)
        logging.info("Searching for %s numbers with base up to %d and prime factor %s.\n", _options.Divides.data(), limit, input.display_text().data());
    else
        logging.info("Searching for %s numbers with prime factor %s.\n", _options.Divides.data(), input.display_text().data());
    if (gwstate.information_only)
        throw TaskAbortException();
    logging.set_prefix("(" + _options.Divides + " div " + input.display_text() + ") ");

    _result.clear();
    Giant Nm1 = *gwstate.N;
    Nm1 -= 1;
    Product Pr(&input, &gwstate, &logging);
    for (int b = 2; b <= limit; b++)
    {
        Base& base = _bases[b - 2];

        auto it = PrimeIterator::get();
        for (; *it < b && b%(*it) != 0; it++);
        if (*it < b)
        {
            Base& base_a = (_bases[*it - 2].power <= _bases[b/(*it) - 2].power ? _bases[*it - 2] : _bases[b/(*it) - 2]);
            Base& base_b = (_bases[*it - 2].power <= _bases[b/(*it) - 2].power ? _bases[b/(*it) - 2] : _bases[*it - 2]);
            base.val = Pr.mul(base_a.val, base_b.val);

            base.sub_val = base_a.sub_val;
            if (base_a.power + _sub <= base_b.power)
                base.sub_val = 1;
            else if (base_a.power < base_b.power)
            {
                for (base.power = base_a.power; base.power < base_b.power; base.power++)
                    base.sub_val = Pr.mul(base.sub_val, base.sub_val);
            }
            base.sub_val = Pr.mul(base.sub_val, base_b.sub_val);
            base.power = base_b.power;
            if (base.sub_val == 1)
            {
                logging.report_param("base_" + base.str, base.power);
                create_task(logging, base);
                if (base.task)
                    std::swap(logging.progress().costs()[logging.progress().cur_stage()], logging.progress().costs().back());
            }
        }
        if (base.sub_val.empty() || base.sub_val == 1)
        {
            run_task(gwstate, file_checkpoint, file_recoverypoint, logging, base);

            _task_check.reset(new CarefulExp(base.exp));
            _task_check->set_error_check(false, true);
            _task_check->init_giant(&input, &gwstate, &logging, std::move(base.sub_val));
            _task_check->run();
            base.sub_val = std::move(_task_check->X0());
            base.val = std::move(*_task_check->result());
        }

        if (base.val == 1 && perfect_power(b) == 1)
        {
            Giant cur_val = base.sub_val;
            int ord = base.power;
            for (; ord < input.factors()[0].second && cur_val != Nm1; ord++)
                cur_val = Pr.mul(cur_val, cur_val);

            if (cur_val != Nm1)
            {
                logging.error("Computation error.\n");
                logging.set_prefix("");
                throw TaskAbortException();
            }

            if (!_result.empty())
                _result += ", ";
            if (b == 2)
                _result += "F(" + std::to_string(ord) + ")";
            else
                _result += "GF(" + std::to_string(ord) + "," + std::to_string(b) + ")";
            if (ord < input.factors()[0].second - 15)
                logging.warning("GF(%d,%d)\n", ord, b);
        }
    }

    if (_options.Divides == "xGF")
    {
        for (int a = 3; a <= limit; a++)
            for (int b = 2; b < a; b++)
                if (gcd(a, b) == 1 && gcd(perfect_power(a), perfect_power(b)) == 1 && (_bases[a - 2].val == _bases[b - 2].val || _bases[a - 2].val + _bases[b - 2].val == *gwstate.N))
                {
                    Base& base_a = (_bases[a - 2].power <= _bases[b - 2].power ? _bases[a - 2] : _bases[b - 2]);
                    Base& base_b = (_bases[a - 2].power <= _bases[b - 2].power ? _bases[b - 2] : _bases[a - 2]);

                    Giant a_val = base_a.sub_val;
                    if (base_a.power + _sub <= base_b.power)
                        a_val = 1;
                    else if (base_a.power < base_b.power)
                    {
                        for (int power = base_a.power; power < base_b.power; power++)
                            a_val = Pr.mul(a_val, a_val);
                    }
                    Giant b_val = base_b.sub_val;

                    int ord = base_b.power;
                    for (; ord < input.factors()[0].second && a_val + b_val != *gwstate.N; ord++)
                    {
                        a_val = Pr.mul(a_val, a_val);
                        b_val = Pr.mul(b_val, b_val);
                    }

                    if (a_val + b_val != *gwstate.N)
                    {
                        Base base(b);
                        base.base.inv(*gwstate.N);
                        base.base *= a;
                        base.base %= (*gwstate.N);
                        base.str = std::to_string(a) + "/" + std::to_string(b);

                        logging.report_param("base_" + base.str, base_b.power);
                        create_task(logging, base);
                        run_task(gwstate, file_checkpoint, file_recoverypoint, logging, base);
                        ord = base.power;
                        b_val = std::move(base.sub_val);

                        for (; ord < input.factors()[0].second && b_val != Nm1; ord++)
                            b_val = Pr.mul(b_val, b_val);

                        if (b_val != Nm1)
                        {
                            logging.error("Computation error.\n");
                            logging.set_prefix("");
                            throw TaskAbortException();
                        }
                    }

                    if (!_result.empty())
                        _result += ", ";
                    _result += "xGF(" + std::to_string(ord) + "," + std::to_string(a) + "," + std::to_string(b) + ")";
                    if (ord < input.factors()[0].second - 15)
                        logging.warning("xGF(%d,%d,%d)\n", ord, a, b);
                }
    }


    logging.set_prefix("");
    if (_result.empty())
    {
        _result = "not found";
        logging.result(false, "%s no divisible numbers found.\n", input.display_text().data());
        logging.result_save(input.input_text() + " no divisible numbers found.\n");
    }
    else
    {
        logging.result(true, "%s divides %s.\n", input.display_text().data(), _result.data());
        logging.result_save(input.input_text() + " divides " + _result + ".\n");
    }

    file_checkpoint.clear(true);
    file_recoverypoint.clear(true);
}
