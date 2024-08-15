#pragma once

#include <set>
#include "fermat.h"

class Pocklington : public Fermat
{
public:
    Pocklington(InputNum& input, Options& options, Logging& logging, Proof* proof);

    void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging, Proof* proof) override;

protected:
    class FactorTask
    {
    public:
        FactorTask(int i) : index(i) { }
        int index;
        std::unique_ptr<CarefulExp> taskFactor;
        std::unique_ptr<CarefulExp> taskCheck;
    };

protected:
    bool _all_factors = false;
    std::vector<FactorTask> _tasks;
    arithmetic::Giant _done;
};

class FactorTree;

class PocklingtonGeneric
{
public:
    PocklingtonGeneric(InputNum& input, Options& options, Logging& logging);

    void run(InputNum& input, arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging);

    bool success() { return _success; }
    bool prime() { return _prime; }
    std::string& res64() { return _res64; }

protected:
    void create_tasks(InputNum& input, Logging& logging, arithmetic::Giant& exp);

protected:
    Options& _options;
    std::unique_ptr<SubLogging> _logging;

    arithmetic::Giant _done;
    std::set<int> _done_factors;
    std::unique_ptr<FactorTree> _tree;
    int _a;

    bool _success = false;
    bool _prime = false;
    std::string _res64;
};

class FactorTree
{
public:
    template<class S>
    FactorTree(S&& exp, int index = -1) : _index(index)
    {
        _exp = std::forward<S>(exp);
    }

    FactorTree(std::vector<std::unique_ptr<FactorTree>>& tree) : _index(-1)
    {
        if (tree.size() == 1)
        {
            _left = std::move(tree[0]);
            return;
        }
        for (auto it = tree.begin(); it != tree.end(); it++)
            (*it)->_left.reset(new FactorTree((*it)->exp(), (*it)->index()));

        while (true)
        {
            std::vector<std::unique_ptr<FactorTree>> next;
            for (auto it = tree.begin(); it != tree.end(); )
            {
                auto& a = *it;  it++;
                if (it == tree.end())
                    next.push_back(std::move(a));
                else
                {
                    auto& b = *it;  it++;
                    swap(a->exp(), b->exp());
                    std::swap(a->_index, b->_index);
                    if (tree.size() == 2)
                    {
                        _left = std::move(tree[0]);
                        _right = std::move(tree[1]);
                        return;
                    }
                    next.emplace_back(new FactorTree(a->exp()*b->exp()));
                    next.back()->_left = std::move(a);
                    next.back()->_right = std::move(b);
                }
            }
            tree = std::move(next);
        }
    }

    int index() { return _index; }
    arithmetic::Giant& exp() { return _exp; }
    bool is_factor() { return !_left && !_right; }
    bool is_last() { return !_right; }
    FactorTree* left() { return _left.get(); }
    FactorTree* right() { return _right.get(); }

private:
    int _index;
    arithmetic::Giant _exp;
    std::unique_ptr<FactorTree> _left;
    std::unique_ptr<FactorTree> _right;
};
