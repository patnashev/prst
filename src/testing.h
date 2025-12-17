#pragma once

int testing_main(int argc, char *argv[]);

struct NTest
{
    uint32_t n;
    uint64_t res64;
    uint64_t cert64;
};

struct BTest
{
    uint32_t b;
    uint64_t res64;
    uint64_t cert64;
};

struct KBNCTest
{
    const char *sk;
    const char *sb;
    uint32_t n;
    int c;
    uint64_t res64;
    uint64_t cert64;
};

struct FreeFormTest
{
    const char *s;
    int bitlen;
    uint64_t res64;
    uint64_t cert64;
};

class Test
{
public:
    Test(int k, BTest& t, int n, int c)
    {
        input.init(k, t.b, n, c);
        input_text = input.input_text();
        input_bitlen = input.bitlen();
        res64 = t.res64;
        cert64 = t.cert64;
    }
    Test(int k, int b, NTest& t, int c)
    {
        input.init(k, b, t.n, c);
        input_text = input.input_text();
        input_bitlen = input.bitlen();
        res64 = t.res64;
        cert64 = t.cert64;
    }
    Test(KBNCTest& t)
    {
        input.init(t.sk, t.sb, t.n, t.c);
        input_text = input.input_text();
        input_bitlen = input.bitlen();
        res64 = t.res64;
        cert64 = t.cert64;
    }
    Test(FreeFormTest& t)
    {
        input_text = t.s;
        input_bitlen = t.bitlen;
        res64 = t.res64 ? t.res64 : 1;
        cert64 = t.cert64;
    }
    virtual ~Test() { }

    double cost() { return input_bitlen*std::sqrt(input_bitlen); }
    virtual std::string display_text()
    {
        if (input.type() == InputNum::ZERO && !input.parse(input_text))
            throw std::runtime_error("Parse failed.");
        if (input.type() == InputNum::KBNC && input.c() == 1 && input.d() == 1 && input.b() == 2 && log2(input.gk()) < input.n())
            return input.display_text() + ", Proth test with certification.";
        return input.display_text() + ", Fermat test with certification.";
    }
    
    virtual void run(Logging& logging, Options& global_options, arithmetic::GWState& global_state);

public:
    std::string input_text;
    int input_bitlen;
    InputNum input;
    uint64_t res64;
    uint64_t cert64;
};

class DeterministicTest : public Test
{
public:
    DeterministicTest(FreeFormTest& t) : Test(t) { }
    DeterministicTest(KBNCTest& t) : Test(t) { }

    std::string display_text() override
    {
        if (input.type() == InputNum::ZERO && !input.parse(input_text))
            throw std::runtime_error("Parse failed.");
        if (input.c() == 1 && (input.type() == InputNum::FACTORIAL || input.type() == InputNum::PRIMORIAL || (input.type() == InputNum::KBNC && (input.n() < 10 || input.factors().size() > 10))))
            return input.display_text() + ", generic Pocklington test.";
        else if (input.c() == 1)
            return input.display_text() + ", Pocklington test.";
        if (input.c() == -1 && (input.type() == InputNum::FACTORIAL || input.type() == InputNum::PRIMORIAL || (input.type() == InputNum::KBNC && (input.n() < 10 || input.factors().size() > 10))))
            return input.display_text() + ", generic Morrison test.";
        else if (input.c() == -1 && (input.b() != 2 || log2(input.gk()) >= input.n()))
            return input.display_text() + ", Morrison test.";
        else if (input.c() == -1)
            return input.display_text() + ", Morrison (LLR) test.";
        return input.display_text() + ", No deterministic test.";
    }

    void run(Logging& logging, Options& global_options, arithmetic::GWState& global_state) override;
};

class TestLogging : public Logging
{
public:
    TestLogging(int level) : Logging(level) { }

    void heartbeat() override
    {
        if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - _last_progress).count() >= Task::PROGRESS_TIME)
        {
            report_progress();
            _last_progress = std::chrono::system_clock::now();
        }
    }

private:
    std::chrono::system_clock::time_point _last_progress = std::chrono::system_clock::now();
};

void RootsTest(Logging& logging, Options& global_options, arithmetic::GWState& global_state);
