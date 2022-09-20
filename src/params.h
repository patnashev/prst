#pragma once

#include <optional>

class Params
{
public:
    std::optional<bool> Check;
    std::optional<bool> CheckNear;

    std::optional<bool> CheckGerbicz;
    std::optional<int> GerbiczCount;
    std::optional<int> GerbiczL;
    std::optional<int> GerbiczL2;

    std::optional<int> SlidingWindow;

    std::optional<int> FermatBase;
    int maxmulbyconst = 1;
};
