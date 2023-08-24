
#include <algorithm>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "order.h"
#include "integer.h"

using namespace arithmetic;

Order::Order(int a, InputNum& input, Params& params, Logging& logging) : _a(a)
{
    _factors = input.factors();
    create_tasks(params, logging, false);

    if (_task)
        params.maxmulbyconst = a;
}

void Order::create_tasks(Params& params, Logging& logging, bool restart)
{
    bool CheckStrong = params.CheckStrong ? params.CheckStrong.value() : false;
    int checks = params.StrongCount ? params.StrongCount.value() : 16;

    _tasks_smooth.clear();
    Giant smooth_base;
    Giant exp;
    exp = 1;
    _task_exp_str.clear();
    for (auto& factor : _factors)
    {
        if (factor.second <= _sub)
            continue;
        if (factor.first == 2 || restart)
        {
            smooth_base = factor.first;
            if (CheckStrong)
                _tasks_smooth.emplace_back(new GerbiczCheckExp(factor.first, factor.second - _sub, checks));
            else
                _tasks_smooth.emplace_back(new SmoothExp(factor.first, factor.second - _sub));
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
        if (CheckStrong && exp.bitlen() > 100)
            _task.reset(new FastLiCheckExp(exp, checks));
        else
            _task.reset(new FastExp(exp));
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

void Order::run(Params& params, InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging)
{
    logging.info("Computing multiplicative order of %d modulo prime %s.\n", _a, input.display_text().data());
    if (gwstate.information_only)
        exit(0);
    logging.set_prefix("ord(" + std::to_string(_a) + ") mod " + input.display_text() + " ");

    while (!_factors.empty())
    {
        Giant sub_val;
        sub_val = _a;

        if (_task)
        {
            logging.info("raising to power %s.\n", _task_exp_str.data());
            if (FastExp* task = dynamic_cast<FastExp*>(_task.get()))
                task->init_small(&input, &gwstate, nullptr, &logging, _a);
            if (FastLiCheckExp* task = dynamic_cast<FastLiCheckExp*>(_task.get()))
                task->init(&input, &gwstate, nullptr, nullptr, &logging, _a);

            _task->run();
            sub_val = *_task->result();
            logging.progress().next_stage();
        }
        for (auto& task_smooth : _tasks_smooth)
        {
            logging.info("raising to power %s^%d.\n", task_smooth->b().to_string().data(), task_smooth->points().back());
            if (SmoothExp* task = dynamic_cast<SmoothExp*>(task_smooth.get()))
                task->init(&input, &gwstate, nullptr, &logging);
            if (GerbiczCheckExp* task = dynamic_cast<GerbiczCheckExp*>(task_smooth.get()))
                task->init(&input, &gwstate, nullptr, nullptr, &logging);
            if (task_smooth->state() == nullptr)
                task_smooth->init_state(new BaseExp::StateValue(0, sub_val));

            task_smooth->run();
            sub_val = *task_smooth->result();
            logging.progress().next_stage();
        }

        if (sub_val == 1)
        {
            for (auto it = _factors.begin(); it != _factors.end(); it++)
                if (it->second <= _sub)
                    it = _factors.erase(it);
                else
                    it->second -= _sub;
            create_tasks(params, logging, true);
            continue;
        }

        logging.info("computing order for each factor.\n");

        _task_check->init(&input, &gwstate, &logging, std::move(sub_val));
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

            factor.task_sub->init(&input, &gwstate, &logging, std::move(sub_val));
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
                factor.task_factor->init(&input, &gwstate, &logging, std::move(cur_val));
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
            create_tasks(params, logging, true);
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
    order += " = (N-1)/";
    order += order_div.to_string();

    logging.set_prefix("");
    logging.result(true, "ord(%d) mod %s = %s.\n", _a, input.display_text().data(), order.data());
    logging.result_save("ord(" + std::to_string(_a) + ") mod " + input.input_text() + " = " + order + ".\n");

    file_checkpoint.clear(true);
    file_recoverypoint.clear(true);
}
