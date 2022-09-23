#pragma once

#include <optional>
#include <string>

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

    std::string ProofPointFilename;
    std::string ProofProductFilename;
    std::optional<int> ProofPointsPerCheck;
    std::optional<int> ProofChecksPerPoint;
    std::string ProofSecuritySeed;
};
