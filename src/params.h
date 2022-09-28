#pragma once

#include <optional>
#include <string>

class Params
{
public:
    int maxmulbyconst = 1;
    int thread_count = 1;
    int spin_threads = 1;

    std::optional<bool> Check;
    std::optional<bool> CheckNear;

    std::optional<bool> CheckGerbicz;
    std::optional<int> GerbiczCount;
    std::optional<int> GerbiczL;
    std::optional<int> GerbiczL2;

    std::optional<int> SlidingWindow;

    std::optional<int> FermatBase;

    std::string ProofPointFilename;
    std::string ProofProductFilename;
    std::optional<int> ProofPointsPerCheck;
    std::optional<int> ProofChecksPerPoint;
    std::string ProofSecuritySeed;

    std::optional<bool> RootOfUnityCheck;
    std::optional<int> RootOfUnitySecurity;
};
