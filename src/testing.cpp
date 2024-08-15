
#include <list>
#include <deque>
#include <tuple>
#include <cmath>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "arithmetic.h"
#include "exception.h"
#include "config.h"
#include "inputnum.h"
#include "file.h"
#include "logging.h"
#include "task.h"
#include "options.h"
#include "fermat.h"
#include "proof.h"
#include "pocklington.h"
#include "morrison.h"
#include "testing.h"

#include "test.data"

using namespace arithmetic;

int testing_main(int argc, char *argv[])
{
    GWState gwstate;
    Options options;
    std::string subset;
    int log_level = Logging::LEVEL_PROGRESS;
    std::string log_file;
    Task::PROGRESS_TIME = 60;

    Config cnfg;
    cnfg.ignore("-test")
        .value_number("-t", 0, gwstate.thread_count, 1, 256)
        .value_number("-t", ' ', gwstate.thread_count, 1, 256)
        .value_number("-spin", ' ', gwstate.spin_threads, 0, 256)
        .value_enum("-cpu", ' ', gwstate.instructions, Enum<std::string>().add("SSE2", "SSE2").add("AVX", "AVX").add("FMA3", "FMA3").add("AVX512F", "AVX512F"))
        .value_number("-fft", '+', gwstate.next_fft_count, 0, 5)
        .group("-fft")
            .value_number("+", 0, gwstate.next_fft_count, 0, 5)
            .value_number("safety", ' ', gwstate.safety_margin, -10.0, 10.0)
            .check("generic", gwstate.force_mod_type, 1)
            .end()
        .group("-check")
            .exclusive()
                .ex_case().check_code("near", [&] { options.CheckNear = true; options.Check = false; }).end()
                .ex_case().check_code("always", [&] { options.CheckNear = false; options.Check = true; }).end()
                .ex_case().check_code("never", [&] { options.CheckNear = false; options.Check = false; }).end()
                .end()
            .group("strong")
                .check("disable", options.CheckStrong, false)
                .value_number("count", ' ', options.StrongCount, 1, 1048576)
                .value_number("L", ' ', options.StrongL, 1, INT_MAX)
                .value_number("L2", ' ', options.StrongL2, 1, INT_MAX)
                .end()
                //.on_check(options.CheckStrong, true)
            .end()
        .group("-factors")
            .check("all", options.AllFactors, true)
            .end()
        .group("-time")
            .value_number("write", ' ', Task::DISK_WRITE_TIME, 1, INT_MAX)
            .value_number("progress", ' ', Task::PROGRESS_TIME, 1, INT_MAX)
            .check("coarse", Task::MULS_PER_STATE_UPDATE, Task::MULS_PER_STATE_UPDATE/10)
            .end()
        .group("-log")
            .exclusive()
                .ex_case().check("debug", log_level, Logging::LEVEL_DEBUG).end()
                .ex_case().check("info", log_level, Logging::LEVEL_INFO).end()
                .ex_case().check("warning", log_level, Logging::LEVEL_WARNING).end()
                .ex_case().check("error", log_level, Logging::LEVEL_ERROR).end()
                .end()
            .value_string("file", ' ', log_file)
            .end()
        .check("-d", log_level, Logging::LEVEL_INFO)
        .value_code("-ini", ' ', [&](const char* param) {
                File ini_file(param, 0);
                ini_file.read_buffer();
                if (ini_file.buffer().empty())
                    printf("ini file not found: %s.\n", param);
                else
                    cnfg.parse_ini(ini_file);
                return true;
            })
        .default_code([&](const char* param) {
                subset = param;
            })
        .parse_args(argc, argv);

    if (subset.empty())
    {
        printf("Usage: PRST -test <subset> <options>\n");
        printf("Subsets:\n");
        printf("\tall = 321plus + 321minus + b5plus + b5minus + gfn13 + special + error + freeform + deterministic + prime\n");
        printf("\tslow = gfn13more + 100186b5minus + 109208b5plus\n");
        return 0;
    }

    TestLogging logging(log_level);
    if (!log_file.empty())
        logging.file_log(log_file);

    std::list<std::tuple<std::string,SubLogging,std::deque<std::unique_ptr<Test>>>> tests;
    auto add = [&](const std::string& subset) -> std::deque<std::unique_ptr<Test>>& { return std::get<2>(tests.emplace_back(subset, SubLogging(logging, log_level), std::deque<std::unique_ptr<Test>>())); };

    if (subset == "all" || subset == "321plus")
    {
        auto& cont = add("321plus");
        for (NTest* nTest = Test321Plus; nTest->n != 0; nTest++)
            cont.emplace_back(new Test(3, 2, *nTest, 1));
    }

    if (subset == "all" || subset == "321minus")
    {
        auto& cont = add("321minus");
        for (NTest* nTest = Test321Minus; nTest->n != 0; nTest++)
            cont.emplace_back(new Test(3, 2, *nTest, -1));
    }

    if (subset == "all" || subset == "b5plus")
    {
        auto& cont = add("b5plus");
        for (NTest* nTest = TestBase5Plus; nTest->n != 0; nTest++)
            cont.emplace_back(new Test(2, 5, *nTest, 1));
    }

    if (subset == "all" || subset == "b5minus")
    {
        auto& cont = add("b5minus");
        for (NTest* nTest = TestBase5Minus; nTest->n != 0; nTest++)
            cont.emplace_back(new Test(2, 5, *nTest, -1));
    }

    if (subset == "all" || subset == "gfn13")
    {
        auto& cont = add("gfn13");
        for (BTest* bTest = TestGFN13; bTest->b != 0; bTest++)
            cont.emplace_back(new Test(1, *bTest, 8192, 1));
    }

    if (subset == "all" || subset == "special")
    {
        auto& cont = add("special");
        for (KBNCTest* kbncTest = TestSpecial; kbncTest->n != 0; kbncTest++)
            cont.emplace_back(new Test(*kbncTest));
    }

    if (subset == "all" || subset == "error")
    {
        add("error");
    }
    
    if (subset == "all" || subset == "freeform")
    {
        auto& cont = add("freeform");
        for (FreeFormTest* ffTest = TestFreeForm; ffTest->s[0] != 0; ffTest++)
        {
            cont.emplace_back(new Test(*ffTest));
            if (ffTest->res64 == 1)
                cont.emplace_back(new DeterministicTest(*ffTest));
        }
    }

    if (subset == "all" || subset == "deterministic")
    {
        std::string str2("2");
        auto& cont = add("deterministic");
        for (KBNCTest* kbncTest = TestPrime; kbncTest->n != 0; kbncTest++)
            if (kbncTest->c == -1 || (kbncTest->c == 1 && kbncTest->sb != str2))
                cont.emplace_back(new DeterministicTest(*kbncTest));
    }

    if (subset == "all" || subset == "prime")
    {
        auto& cont = add("prime");
        for (KBNCTest* kbncTest = TestPrime; kbncTest->n != 0; kbncTest++)
            cont.emplace_back(new Test(*kbncTest));
    }
    
    if (subset == "slow" || subset == "gfn13more")
    {
        auto& cont = add("gfn13more");
        for (BTest* bTest = TestGFN13More; bTest->b != 0; bTest++)
            cont.emplace_back(new Test(1, *bTest, 8192, 1));
    }

    if (subset == "slow" || subset == "109208b5plus")
    {
        auto& cont = add("109208b5plus");
        for (NTest* nTest = Test109208Base5Plus; nTest->n != 0; nTest++)
            cont.emplace_back(new Test(109208, 5, *nTest, 1));
    }

    if (subset == "slow" || subset == "100186b5minus")
    {
        auto& cont = add("100186b5minus");
        for (NTest* nTest = Test100186Base5Minus; nTest->n != 0; nTest++)
            cont.emplace_back(new Test(100186, 5, *nTest, -1));
    }

    for (auto& subsetTests : tests)
    {
        for (auto& test : std::get<2>(subsetTests))
            std::get<1>(subsetTests).progress().add_stage(test->cost());
        logging.progress().add_stage(std::get<1>(subsetTests).progress().cost_total());
    }

    const char* test_text = "";
    try
    {
        for (auto& subsetTests : tests)
        {
            logging.progress().update(0, 0);
            logging.warning("Running %s tests.\n", std::get<0>(subsetTests).data());
            if (std::get<0>(subsetTests) == "error")
            {
                test_text = "error";
                SubLogging subLogging(std::get<1>(subsetTests), log_level > Logging::LEVEL_DEBUG ? Logging::LEVEL_ERROR : Logging::LEVEL_INFO);
                if (log_level > Logging::LEVEL_DEBUG)
                    subLogging.file_result(log_file);
                else
                    subLogging.file_log(log_file);
                RootsTest(subLogging, options, gwstate);
            }
            else
                for (auto& test : std::get<2>(subsetTests))
                {
                    logging.info("%s\n", test->display_text().data());
                    std::get<1>(subsetTests).progress().update(0, 0);
                    SubLogging subLogging(std::get<1>(subsetTests), log_level > Logging::LEVEL_DEBUG ? Logging::LEVEL_ERROR : Logging::LEVEL_INFO);
                    if (log_level > Logging::LEVEL_DEBUG)
                        subLogging.file_result(log_file);
                    else
                        subLogging.file_log(log_file);
                    test->run(subLogging, options, gwstate);
                    std::get<1>(subsetTests).progress().next_stage();
                }
            logging.progress().next_stage();
        }
        logging.error("All tests completed successfully.\n");
    }
    catch (const TaskAbortException&)
    {
        if (Task::abort_flag())
            logging.warning("Test aborted.\n");
        else
            logging.error("Test %s failed.\n", test_text);
    }

    return 0;
}

void Test::run(Logging& logging, Options& global_options, GWState& global_state)
{
    Options options;
    options.Check = global_options.Check;
    options.CheckNear = global_options.CheckNear;
    options.CheckStrong = global_options.CheckStrong;
    options.ProofSecuritySeed = "12345";
    options.RootOfUnityCheck = false;
    int proof_count = 16;

    uint32_t fingerprint = input.fingerprint();
    File file_cert("prst_cert", fingerprint);
    Proof proof(Proof::SAVE, proof_count, input, options, file_cert, logging);
    Fermat fermat(Fermat::AUTO, input, options, logging, &proof);

    fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat.a()) + "." + std::to_string(proof.points()[proof_count].pos));
    File file_proofpoint("prst_proof", fingerprint);
    File file_proofproduct("prst_prod", fingerprint);
    File file_checkpoint("prst_ckpt", fingerprint);
    File file_recoverypoint("prst_rcpt", fingerprint);

    GWState gwstate;
    gwstate.copy(global_state);
    gwstate.maxmulbyconst = options.maxmulbyconst;
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());

    auto finally = [&]
    {
        file_cert.clear();
        file_proofpoint.clear(true);
        file_proofproduct.clear(true);
        file_checkpoint.clear(true);
        file_recoverypoint.clear(true);
        gwstate.done();
    };
    try
    {
        proof.init_files(&file_proofpoint, &file_proofproduct, &file_cert);
        fermat.run(input, gwstate, file_checkpoint, file_recoverypoint, logging, &proof);
        if (res64 == 0)
            res64 = fermat.success() ? 1 : std::stoull(fermat.res64(), nullptr, 16);
        if (cert64 == 0)
            cert64 = std::stoull(proof.res64(), nullptr, 16);
        if (fermat.success() != (res64 == 1))
        {
            logging.error("Primality mismatch.\n");
            throw TaskAbortException();
        }
        if (!fermat.success() && std::stoull(fermat.res64(), nullptr, 16) != res64)
        {
            logging.error("RES64 mismatch.\n");
            throw TaskAbortException();
        }
        if (std::stoull(proof.res64(), nullptr, 16) != cert64)
        {
            logging.error("Raw certificate mismatch.\n");
            throw TaskAbortException();
        }

        gwstate.done();
        gwstate.next_fft_count = 1;
        input.setup(gwstate);
        logging.info("Using %s.\n", gwstate.fft_description.data());

        Proof proof_build(Proof::BUILD, proof_count, input, options, file_cert, logging);
        proof_build.calc_points(proof.points()[proof_count].pos, !proof.Li(), input, options, logging);
        proof_build.init_files(&file_proofpoint, &file_proofproduct, &file_cert);
        logging.progress().add_stage(1);
        logging.progress().add_stage(proof_build.cost());
        fermat.run(input, gwstate, file_checkpoint, file_recoverypoint, logging, &proof_build);
        if (std::stoull(dynamic_cast<ProofBuild*>(proof_build.task())->raw_res64(), nullptr, 16) != cert64)
        {
            logging.error("Build raw certificate mismatch.\n");
            throw TaskAbortException();
        }
        if (fermat.success() != (res64 == 1))
        {
            logging.error("Build primality mismatch.\n");
            throw TaskAbortException();
        }
        if (!fermat.success() && std::stoull(fermat.res64(), nullptr, 16) != res64)
        {
            logging.error("Build RES64 mismatch.\n");
            throw TaskAbortException();
        }

        gwstate.done();
        gwstate.next_fft_count = 0;
        input.setup(gwstate);
        logging.info("Using %s.\n", gwstate.fft_description.data());

        Proof proof_cert(Proof::CERT, 0, input, options, file_cert, logging);
        proof_cert.run(input, gwstate, file_checkpoint, file_recoverypoint, logging);
        if (proof_cert.res64() != proof_build.res64())
        {
            logging.error("Certificate mismatch.\n");
            throw TaskAbortException();
        }
    }
    catch (const TaskAbortException&)
    {
        finally();
        throw;
    }
    finally();
}

void DeterministicTest::run(Logging& logging, Options& global_options, GWState& global_state)
{
    Options options;
    options.Check = global_options.Check;
    options.CheckNear = global_options.CheckNear;
    options.CheckStrong = global_options.CheckStrong;
    options.AllFactors = global_options.AllFactors;

    uint32_t fingerprint = input.fingerprint();
    input.factorize_f_p();
    if (!input.is_half_factored())
        throw std::runtime_error("Not enough factors.");

    std::unique_ptr<Pocklington> pocklington;
    std::unique_ptr<PocklingtonGeneric> pocklingtonGeneric;
    std::unique_ptr<Morrison> morrison;
    if (input.c() == 1 && (input.type() == InputNum::FACTORIAL || input.type() == InputNum::PRIMORIAL || (input.type() == InputNum::KBNC && input.bitlen() > 1000 && input.n() < 10)))
        pocklingtonGeneric.reset(new PocklingtonGeneric(input, options, logging));
    else if (input.c() == 1)
        pocklington.reset(new Pocklington(input, options, logging, nullptr));
    if (input.c() == -1 && (input.type() == InputNum::FACTORIAL || input.type() == InputNum::PRIMORIAL || (input.type() == InputNum::KBNC && input.bitlen() > 1000 && input.n() < 10)))
        morrison.reset(new MorrisonGeneric(input, options, logging));
    else if (input.c() == -1)
        morrison.reset(new Morrison(input, options, logging));

    File file_checkpoint("prst_ckpt", fingerprint);
    File file_recoverypoint("prst_rcpt", fingerprint);
    File file_params("prst_param", fingerprint);
    logging.file_progress(&file_params);

    GWState gwstate;
    gwstate.copy(global_state);
    gwstate.maxmulbyconst = options.maxmulbyconst;
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());

    auto finally = [&]
    {
        file_checkpoint.clear(true);
        file_recoverypoint.clear(true);
        file_params.clear(true);
        gwstate.done();
    };
    try
    {
        bool prime = false;
        std::string sres64 = "0";
        if (pocklington)
        {
            pocklington->run(input, gwstate, file_checkpoint, file_recoverypoint, logging, nullptr);
            prime = pocklington->prime();
            sres64 = pocklington->res64();
        }
        if (pocklingtonGeneric)
        {
            pocklingtonGeneric->run(input, gwstate, file_checkpoint, file_recoverypoint, logging);
            prime = pocklingtonGeneric->prime();
            sres64 = pocklingtonGeneric->res64();
        }
        if (morrison)
        {
            morrison->run(input, gwstate, file_checkpoint, file_recoverypoint, logging);
            prime = morrison->prime();
            sres64 = morrison->res64();
        }
        if (res64 == 0)
            res64 = prime ? 1 : std::stoull(sres64, nullptr, 16);
        if (prime != (res64 == 1))
        {
            logging.error("Primality mismatch.\n");
            throw TaskAbortException();
        }
        if (!prime && std::stoull(sres64, nullptr, 16) != res64)
        {
            logging.error("RES64 mismatch.\n");
            throw TaskAbortException();
        }
    }
    catch (const TaskAbortException&)
    {
        finally();
        throw;
    }
    finally();
}

void RootsTest(Logging& logging, Options& global_options, GWState& global_state)
{
    SubLogging noLogging(logging, Logging::LEVEL_ERROR + 1);
    Options options;
    options.Check = global_options.Check;
    options.CheckNear = global_options.CheckNear;
    options.CheckStrong = global_options.CheckStrong;
    options.RootOfUnityCheck = false;
    InputNum input;
    int proof_count = 4;

    input.parse("3*2^353+1");

    uint32_t fingerprint = input.fingerprint();
    File file_cert("prst_cert", fingerprint);
    Proof proof(Proof::SAVE, proof_count, input, options, file_cert, logging);
    Fermat fermat(Fermat::AUTO, input, options, logging, &proof);
    Proof proof_build(Proof::BUILD, proof_count, input, options, file_cert, logging);
    proof_build.calc_points(proof.points()[proof_count].pos, !proof.Li(), input, options, logging);

    fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat.a()) + "." + std::to_string(proof.points()[proof_count].pos));
    File file_proofpoint("prst_proof", fingerprint);
    File file_proofproduct("prst_prod", fingerprint);
    File file_checkpoint("prst_ckpt", fingerprint);
    File file_recoverypoint("prst_rcpt", fingerprint);

    proof.init_files(&file_proofpoint, &file_proofproduct, &file_cert);
    proof_build.init_files(&file_proofpoint, &file_proofproduct, &file_cert);

    GWState gwstate;
    gwstate.copy(global_state);
    gwstate.maxmulbyconst = fermat.a();
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());
    GWArithmetic gw(gwstate);

    auto finally = [&]
    {
        file_cert.clear();
        file_proofpoint.clear(true);
        file_proofproduct.clear(true);
        file_checkpoint.clear(true);
        file_recoverypoint.clear(true);
        gwstate.done();
    };

    std::vector<std::unique_ptr<BaseExp::State>> points;
    points.emplace_back();
    try
    {
        fermat.run(input, gwstate, file_checkpoint, file_recoverypoint, logging, &proof);

        for (int i = 1; i <= proof_count; i++)
            points.emplace_back(BaseExp::State::read_file(file_proofpoint.children()[i].get()));
        
        GWNum X(gw);
        X = Giant::rnd(gwstate.N->bitlen());
        BaseExp::StateSerialized point;
        point.set(points[1]->iteration(), X);
        file_proofpoint.children()[1]->write(point);

        file_proofproduct.clear(true);
        proof.run(input, gwstate, logging, nullptr);
        fermat.run(input, gwstate, file_checkpoint, file_recoverypoint, logging, &proof_build);

        Proof proof_cert(Proof::CERT, 0, input, options, file_cert, logging);
        proof_cert.run(input, gwstate, file_checkpoint, file_recoverypoint, logging);
        if (proof_cert.res64() == proof_build.res64())
        {
            logging.error("Certificate failure.\n");
            throw TaskAbortException();
        }
    }
    catch (const TaskAbortException&)
    {
        finally();
        throw;
    }

    auto test_attack = [&](Fermat& fermat)
    {
        proof.run(input, gwstate, logging, nullptr);
        fermat.run(input, gwstate, file_checkpoint, file_recoverypoint, logging, &proof_build);
        if (fermat.success())
        {
            logging.error("Attack failure.\n");
            throw TaskAbortException();
        }

        Proof proof_cert(Proof::CERT, 0, input, options, file_cert, logging);
        proof_cert.run(input, gwstate, file_checkpoint, file_recoverypoint, logging);
        if (proof_cert.res64() != proof_build.res64())
        {
            logging.error("Attacking certificate mismatch.\n");
            throw TaskAbortException();
        }

        Proof proof_root(Proof::ROOT, proof_count, input, options, file_cert, logging);
        try
        {
            proof_root.run(input, gwstate, logging.level() > Logging::LEVEL_INFO ? noLogging : logging, &fermat.result());
        }
        catch (const TaskAbortException&) {}
        if (*proof_root.taskRoot()->result() != 1)
        {
            logging.error("Roots of unity check failed to detect the attack.\n");
            throw TaskAbortException();
        }
    };

    try
    {
        Giant tmp;
        tmp = 1;
        FastExp task(tmp << input.n());
        task.init(&input, &gwstate, nullptr, &logging, fermat.a());
        task.run();
        GWNum R(gw);
        R = *task.result();

        GWNum X(gw);
        points[4]->to_GWNum(X);
        X *= R;
        BaseExp::StateValue pointv;
        pointv.set(points[4]->iteration(), X);
        file_proofpoint.children()[4]->write(pointv);

        R.square();
        points[1]->to_GWNum(X);
        X *= R;
        BaseExp::StateSerialized point;
        point.set(points[1]->iteration(), X);
        file_proofpoint.children()[1]->write(point);

        file_proofproduct.clear(true);
        test_attack(fermat);
    }
    catch (const TaskAbortException&)
    {
        finally();
        throw;
    }

    try
    {
        FastExp task(input.gk() << (input.n() - 2));
        task.init(&input, &gwstate, nullptr, &logging, fermat.a());
        task.run();
        GWNum R(gw);
        R = *task.result();

        GWNum X(gw);
        points[4]->to_GWNum(X);
        X *= R;
        BaseExp::StateValue pointv;
        pointv.set(points[4]->iteration(), X);
        file_proofpoint.children()[4]->write(pointv);

        points[2]->to_GWNum(X);
        X *= R;
        BaseExp::StateSerialized point;
        point.set(points[2]->iteration(), X);
        file_proofpoint.children()[2]->write(point);

        file_proofpoint.children()[1]->write(*points[1]);

        file_proofproduct.clear(true);
        test_attack(fermat);
    }
    catch (const TaskAbortException&)
    {
        finally();
        throw;
    }
    finally();

    input.parse("960^128+1");

    proof = Proof(Proof::SAVE, proof_count, input, options, file_cert, logging);
    Fermat fermat2(Fermat::AUTO, input, options, logging, &proof);
    proof_build = Proof(Proof::BUILD, proof_count, input, options, file_cert, logging);
    proof_build.calc_points(proof.points()[proof_count].pos, !proof.Li(), input, options, logging);

    fingerprint = input.fingerprint();
    file_cert = File("prst_cert", fingerprint);
    fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat2.a()) + "." + std::to_string(proof.points()[proof_count].pos));
    file_proofpoint = File("prst_proof", fingerprint);
    file_proofproduct = File("prst_prod", fingerprint);
    file_checkpoint = File("prst_ckpt", fingerprint);
    file_recoverypoint = File("prst_rcpt", fingerprint);

    proof.init_files(&file_proofpoint, &file_proofproduct, &file_cert);
    proof_build.init_files(&file_proofpoint, &file_proofproduct, &file_cert);

    gwstate.copy(global_state);
    gwstate.maxmulbyconst = fermat2.a();
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());

    try
    {
        fermat2.run(input, gwstate, file_checkpoint, file_recoverypoint, logging, &proof);

        FastExp task(power(input.gb()/3, 3)*power(input.gb(), input.n() - 3));
        task.init(&input, &gwstate, nullptr, &logging, fermat2.a());
        task.run();
        GWNum R(gw);
        R = *task.result();
        R.square();

        BaseExp::StateValue pointv;
        file_proofpoint.children()[4]->read(pointv);
        GWNum X(gw);
        pointv.to_GWNum(X);
        X *= R;
        pointv.set(pointv.iteration(), X);
        file_proofpoint.children()[4]->write(pointv);

        std::unique_ptr<BaseExp::State> point;
        point.reset(BaseExp::State::read_file(file_proofpoint.children()[2].get()));
        point->to_GWNum(X);
        X *= R;
        point->set(point->iteration(), X);
        file_proofpoint.children()[2]->write(*point);

        R.square();
        point.reset(BaseExp::State::read_file(file_proofpoint.children()[3].get()));
        point->to_GWNum(X);
        X *= R;
        X *= R;
        X *= R;
        X *= R;
        X *= R;
        point->set(point->iteration(), X);
        file_proofpoint.children()[3]->write(*point);

        file_proofproduct.clear(true);
        test_attack(fermat2);
    }
    catch (const TaskAbortException&)
    {
        finally();
        throw;
    }
    finally();

    input.parse("2*5^178-1");

    proof = Proof(Proof::SAVE, proof_count, input, options, file_cert, logging);
    Fermat fermat3(Fermat::AUTO, input, options, logging, &proof);
    proof_build = Proof(Proof::BUILD, proof_count, input, options, file_cert, logging);
    proof_build.calc_points(proof.points()[proof_count].pos, !proof.Li(), input, options, logging);

    fingerprint = input.fingerprint();
    file_cert = File("prst_cert", fingerprint);
    fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat3.a()) + "." + std::to_string(proof.points()[proof_count].pos));
    file_proofpoint = File("prst_proof", fingerprint);
    file_proofproduct = File("prst_prod", fingerprint);
    file_checkpoint = File("prst_ckpt", fingerprint);
    file_recoverypoint = File("prst_rcpt", fingerprint);

    proof.init_files(&file_proofpoint, &file_proofproduct, &file_cert);
    proof_build.init_files(&file_proofpoint, &file_proofproduct, &file_cert);

    gwstate.copy(global_state);
    gwstate.maxmulbyconst = fermat3.a();
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());

    try
    {
        fermat3.run(input, gwstate, file_checkpoint, file_recoverypoint, logging, &proof);

        FastExp task((input.value() - 1)/3);
        task.init(&input, &gwstate, nullptr, &logging, 2);
        task.run();
        GWNum R(gw);
        R = *task.result();

        BaseExp::StateValue pointv;
        file_proofpoint.children()[4]->read(pointv);
        GWNum X(gw);
        pointv.to_GWNum(X);
        X *= R;
        pointv.set(pointv.iteration(), X);
        file_proofpoint.children()[4]->write(pointv);

        R.square();
        std::unique_ptr<BaseExp::State> point(BaseExp::State::read_file(file_proofpoint.children()[1].get()));
        point->to_GWNum(X);
        X *= R;
        point->set(point->iteration(), X);
        file_proofpoint.children()[1]->write(*point);

        file_proofproduct.clear(true);
        test_attack(fermat3);
    }
    catch (const TaskAbortException&)
    {
        finally();
        throw;
    }
    finally();
}
