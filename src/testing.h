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

class Test
{
public:
    Test(int k, BTest& t, int n, int c)
    {
        input.init(k, t.b, n, c);
        res64 = t.res64;
        cert64 = t.cert64;
    }
    Test(int k, int b, NTest& t, int c)
    {
        input.init(k, b, t.n, c);
        res64 = t.res64;
        cert64 = t.cert64;
    }
    Test(KBNCTest& t)
    {
        input.init(t.sk, t.sb, t.n, t.c);
        res64 = t.res64;
        cert64 = t.cert64;
    }
    virtual ~Test() { }

    double cost() { double len = log2(input.gb())*input.n(); return len*std::sqrt(len)*(input.b() != 2 ? 1.2 : 1.0); }
    virtual void run(Logging& logging, Options& global_options, arithmetic::GWState& global_state);

public:
    InputNum input;
    uint64_t res64;
    uint64_t cert64;
};

class DeterministicTest : public Test
{
public:
    DeterministicTest(KBNCTest& t) : Test(t) { }

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
