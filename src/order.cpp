
#include <algorithm>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "order.h"
#include "integer.h"

using namespace arithmetic;

Order::Order(InputNum& a, InputNum& input, Params& params, Logging& logging)
{
    _factors = input.factors();
    Giant ga = a.value();
    create_tasks(ga, params, logging, false);

    if (_task && ga <= GWMULBYCONST_MAX)
        params.maxmulbyconst = ga.data()[0];
}

void Order::create_tasks(Giant& a, Params& params, Logging& logging, bool restart)
{
    bool CheckStrong = params.CheckStrong ? params.CheckStrong.value() : true;
    int checks = params.StrongCount ? params.StrongCount.value() : 16;

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
            if (CheckStrong)
            {
                _tasks_smooth.emplace_back(new GerbiczCheckExp(factor.first, factor.second - _sub, checks, std::bind(&Order::on_point, this, std::placeholders::_1, std::placeholders::_2)));
                for (auto& p : _tasks_smooth.back()->points())
                    p.value = true;
            }
            else
            {
                int n = factor.second - _sub;
                int len = n/checks;
                if (len == 0)
                    len = 1;
                std::vector<MultipointExp::Point> points;
                for (int i = 0; i <= checks && len*i <= n; i++)
                    points.emplace_back(len*i);
                if (points.back().pos != n)
                    points.emplace_back(n);
                _tasks_smooth.emplace_back(new MultipointExp(factor.first, true, points, std::bind(&Order::on_point, this, std::placeholders::_1, std::placeholders::_2)));
            }
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
}

void Order::run(InputNum& a, Params& params, InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging)
{
    Giant ga = a.value();
    logging.info("Computing multiplicative order of %s modulo prime %s.\n", a.display_text().data(), input.display_text().data());
    if (gwstate.information_only)
        exit(0);
    logging.set_prefix("ord(" + a.display_text() + ") mod " + input.display_text() + " ");

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
                it->second = task_smooth->points()[_task_break].pos;
            }
            logging.progress().next_stage();
        }

        if (sub_val == 1)
        {
            for (auto it = _factors.begin(); it != _factors.end(); it++)
                if (it->second > _sub)
                    it->second -= _sub;
            create_tasks(ga, params, logging, true);
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
            create_tasks(ga, params, logging, true);
    }

    Giant order_div;
    order_div = 1;
    std::string order;
    for (auto& factor : _order)
    {
        auto it = input.factors().begin();
        for (; it != input.factors().end() && it->first != factor.first; it++);
        if (it->second != factor.second)
            order_div *= power(factor.first, it->second - factor.second);

        if (!order.empty())
            order += "*";
        order += factor.first.to_string();
        if (factor.second > 1)
        {
            order += "^";
            order += std::to_string(factor.second);
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
        order += " = (N-1)/";
        order += order_div.to_string();
    }

    logging.set_prefix("");
    logging.result(true, "ord(%s) mod %s = %s.\n", a.display_text().data(), input.display_text().data(), order.data());
    logging.result_save("ord(" + a.display_text() + ") mod " + input.input_text() + " = " + order + ".\n");

    file_checkpoint.clear(true);
    file_recoverypoint.clear(true);
}
