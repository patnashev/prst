#pragma once

#include <optional>
#include <string>

#define PRST_EXIT_NORMAL 0
#define PRST_EXIT_FAILURE 1
#define PRST_EXIT_PRIMEFOUND 2

class Options
{
public:
    
    // Modes

    bool ForceFermat = false;
    std::unique_ptr<InputNum> OrderA;

    // Options

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
    std::optional<bool> ProofPointWriteMode;
    std::string ProofSecuritySeed;

    std::optional<bool> RootOfUnityCheck;
    std::optional<int> RootOfUnitySecurity;

    std::optional<bool> AllFactors;

    // Set by Run::create()
    
    int maxmulbyconst = 1;
};

class Proof;

class Run
{
public:
    Run(InputNum& input_, Options& options) : input(input_), _options(options) { }
    Run(const char* name, InputNum& input_, Options& options) : _name(name), input(input_), _options(options) { }
    virtual ~Run() { }

    const std::string& name() { return _name; }
    uint32_t fingerprint() { return _fingerprint; }
    bool success() { return _success; }
    bool prime() { return _prime; }
    std::string& res64() { return _res64; }

    virtual void run(arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging) = 0;

    static Run* create(InputNum& input, Options& options, Logging& logging, Proof* proof = nullptr);

protected:
    InputNum& input;
    Options& _options;
    std::string _name;
    uint32_t _fingerprint = 0;
    bool _success = false;
    bool _prime = false;
    std::string _res64;
};
