#pragma once

#include <optional>
#include <string>

#define PRST_EXIT_NORMAL 0
#define PRST_EXIT_FAILURE 1
#define PRST_EXIT_PRIMEFOUND 2

class Options
{
public:
    int maxmulbyconst = 1;

    std::optional<bool> Check;
    std::optional<bool> CheckNear;

    std::optional<bool> CheckStrong;
    std::optional<int> StrongCount;
    std::optional<int> StrongL;
    std::optional<int> StrongL2;

    std::optional<int> SlidingWindow;

    std::optional<int> FermatBase;

    std::string ProofPointFilename;
    std::string ProofProductFilename;
    std::optional<int> ProofPointsPerCheck;
    std::optional<int> ProofChecksPerPoint;
    std::string ProofSecuritySeed;

    std::optional<bool> RootOfUnityCheck;
    std::optional<int> RootOfUnitySecurity;

    std::optional<bool> AllFactors;
};
