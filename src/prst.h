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
    std::string Divides;
    std::optional<int> DividesLimit;

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

    // GWState options

    int thread_count = 1;
    int spin_threads = 1;
    std::string instructions;
    int next_fft_count = 0;
    double safety_margin = 0;
    int force_mod_type = 0;
    bool information_only = false;

    void configure(arithmetic::GWState& gwstate)
    {
        gwstate.thread_count = thread_count;
        gwstate.next_fft_count = next_fft_count;
        gwstate.safety_margin = safety_margin;
        gwstate.force_mod_type = force_mod_type;
        gwstate.spin_threads = spin_threads;
        gwstate.instructions = instructions;

        gwstate.information_only = information_only;
        gwstate.known_factors = 1;
    }
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
    arithmetic::Giant& factor() { return _factor; }
    std::string& result() { return _result; }
    bool finished() { return _success || _prime || !_res64.empty() || !_factor.empty() || !_result.empty(); }

    static Run* create(InputNum& input, Options& options, Logging& logging, Proof* proof = nullptr);
    virtual void run(arithmetic::GWState& gwstate, File& file_checkpoint, File& file_recoverypoint, Logging& logging) = 0;
    void primality_result(Logging& logging);
    static void result_prime(InputNum& input, Logging& logging, double time);
    static void result_probable_prime(InputNum& input, Logging& logging, double time);
    static void result_not_prime(InputNum& input, Logging& logging, double time);
    static void result_not_prime_divisible(InputNum& input, Logging& logging, arithmetic::Giant& factor, double time);
    static void result_not_prime_res64(InputNum& input, Logging& logging, const std::string& res64, double time);
    static void result_not_prime_but_probable_res64(InputNum& input, Logging& logging, const std::string& res64, double time);
    static void result_not_probable_prime_res64(InputNum& input, Logging& logging, const std::string& res64, double time);

protected:
    InputNum& input;
    Options& _options;
    std::string _name;
    uint32_t _fingerprint = 0;
    bool _success = false;
    bool _prime = false;
    std::string _res64;
    arithmetic::Giant _factor;
    std::string _result;
};
