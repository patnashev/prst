
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
#include "params.h"
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
    Params params;
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
            .check("generic", gwstate.force_general_mod, true)
            .end()
        .group("-check")
            .exclusive()
                .ex_case().check_code("near", [&] { params.CheckNear = true; params.Check = false; }).end()
                .ex_case().check_code("always", [&] { params.CheckNear = false; params.Check = true; }).end()
                .ex_case().check_code("never", [&] { params.CheckNear = false; params.Check = false; }).end()
                .end()
            .group("strong")
                .value_number("count", ' ', params.StrongCount, 1, 1048576)
                .value_number("L", ' ', params.StrongL, 1, INT_MAX)
                .value_number("L2", ' ', params.StrongL2, 1, INT_MAX)
                .end()
                .on_check(params.CheckStrong, true)
            .end()
        .group("-factors")
            .check("all", params.AllFactors, true)
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
        printf("\tall = 321plus + 321minus + b5plus + b5minus + gfn13 + special + error + deterministic + prime\n");
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
        auto& cont = add("error");
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
                RootsTest(subLogging, params, gwstate);
            }
            else
                for (auto& test : std::get<2>(subsetTests))
                {
                    test_text = test->input.display_text().data();
                    logging.info("%s\n", test_text);
                    std::get<1>(subsetTests).progress().update(0, 0);
                    SubLogging subLogging(std::get<1>(subsetTests), log_level > Logging::LEVEL_DEBUG ? Logging::LEVEL_ERROR : Logging::LEVEL_INFO);
                    if (log_level > Logging::LEVEL_DEBUG)
                        subLogging.file_result(log_file);
                    else
                        subLogging.file_log(log_file);
                    test->run(subLogging, params, gwstate);
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

void Test::run(Logging& logging, Params& global_params, GWState& global_state)
{
    Params params;
    params.Check = global_params.Check;
    params.CheckNear = global_params.CheckNear;
    params.CheckStrong = global_params.CheckStrong;
    params.ProofSecuritySeed = "12345";
    params.RootOfUnityCheck = false;
    int proof_count = 16;

    uint32_t fingerprint = input.fingerprint();
    File file_cert("prst_cert", fingerprint);
    Proof proof(Proof::SAVE, proof_count, input, params, file_cert, logging);
    Fermat fermat(Fermat::AUTO, input, params, logging, &proof);

    fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat.a()) + "." + std::to_string(proof.points()[proof_count]));
    File file_proofpoint("prst_proof", fingerprint);
    File file_proofproduct("prst_prod", fingerprint);
    File file_checkpoint("prst_c", fingerprint);
    File file_recoverypoint("prst_r", fingerprint);

    GWState gwstate;
    gwstate.copy(global_state);
    gwstate.maxmulbyconst = params.maxmulbyconst;
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

        Proof proof_build(Proof::BUILD, proof_count, input, params, file_cert, logging);
        proof_build.points() = std::move(proof.points());
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

        Proof proof_cert(Proof::CERT, 0, input, params, file_cert, logging);
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

void DeterministicTest::run(Logging& logging, Params& global_params, GWState& global_state)
{
    Params params;
    params.Check = global_params.Check;
    params.CheckNear = global_params.CheckNear;
    params.CheckStrong = global_params.CheckStrong;
    params.AllFactors = global_params.AllFactors;

    uint32_t fingerprint = input.fingerprint();
    std::unique_ptr<Pocklington> pocklington;
    std::unique_ptr<Morrison> morrison;
    if (input.c() == 1)
        pocklington.reset(new Pocklington(input, params, logging, nullptr));
    if (input.c() == -1)
        morrison.reset(new Morrison(input, params, logging));
    File file_checkpoint("prst_c", fingerprint);
    File file_recoverypoint("prst_r", fingerprint);
    File file_params("prst_p", fingerprint);

    GWState gwstate;
    gwstate.copy(global_state);
    gwstate.maxmulbyconst = params.maxmulbyconst;
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
        if (morrison)
        {
            morrison->run(input, gwstate, file_checkpoint, file_params, logging);
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

void RootsTest(Logging& logging, Params& global_params, GWState& global_state)
{
    SubLogging noLogging(logging, Logging::LEVEL_ERROR + 1);
    Params params;
    params.Check = global_params.Check;
    params.CheckNear = global_params.CheckNear;
    params.CheckStrong = global_params.CheckStrong;
    params.RootOfUnityCheck = false;
    InputNum input;
    int proof_count = 4;
    Giant tmp, root;
    BaseExp::State point1;
    BaseExp::State point2;
    BaseExp::State point3;
    BaseExp::State point4;

    input.parse("3*2^353+1");

    uint32_t fingerprint = input.fingerprint();
    File file_cert("prst_cert", fingerprint);
    Proof proof(Proof::SAVE, proof_count, input, params, file_cert, logging);
    Fermat fermat(Fermat::AUTO, input, params, logging, &proof);
    Proof proof_build(Proof::BUILD, proof_count, input, params, file_cert, logging);
    proof_build.points() = proof.points();

    fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat.a()) + "." + std::to_string(proof.points()[proof_count]));
    File file_proofpoint("prst_proof", fingerprint);
    File file_proofproduct("prst_prod", fingerprint);
    File file_checkpoint("prst_c", fingerprint);
    File file_recoverypoint("prst_r", fingerprint);

    proof.init_files(&file_proofpoint, &file_proofproduct, &file_cert);
    proof_build.init_files(&file_proofpoint, &file_proofproduct, &file_cert);

    GWState gwstate;
    gwstate.copy(global_state);
    gwstate.maxmulbyconst = fermat.a();
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
        fermat.run(input, gwstate, file_checkpoint, file_recoverypoint, logging, &proof);

        file_proofpoint.children()[1]->read(point1);
        tmp = std::move(point1.X());
        point1.X() = Giant::rnd(gwstate.N->bitlen());
        file_proofpoint.children()[1]->write(point1);
        point1.X() = std::move(tmp);
        file_proofproduct.clear(true);

        proof.run(input, gwstate, logging, nullptr);
        fermat.run(input, gwstate, file_checkpoint, file_recoverypoint, logging, &proof_build);

        Proof proof_cert(Proof::CERT, 0, input, params, file_cert, logging);
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

        Proof proof_cert(Proof::CERT, 0, input, params, file_cert, logging, false);
        proof_cert.run(input, gwstate, file_checkpoint, file_recoverypoint, logging);
        if (proof_cert.res64() != proof_build.res64())
        {
            logging.error("Attacking certificate mismatch.\n");
            throw TaskAbortException();
        }

        Proof proof_root(Proof::ROOT, proof_count, input, params, file_cert, logging);
        try
        {
            proof_root.run(input, gwstate, logging.level() > Logging::LEVEL_INFO ? noLogging : logging, &fermat.result());
        }
        catch (const TaskAbortException&) {}
        if (proof_root.taskRoot()->state()->X() != 1)
        {
            logging.error("Roots of unity check failed to detect the attack.\n");
            throw TaskAbortException();
        }
    };

    try
    {
        tmp = 1;
        FastExp task(tmp << input.n());
        task.init(&input, &gwstate, nullptr, &logging, fermat.a());
        task.run();
        root = std::move(task.state()->X());

        file_proofpoint.children()[4]->read(point4);
        tmp = std::move(point4.X());
        point4.X() = tmp*root%*gwstate.N;
        file_proofpoint.children()[4]->write(point4);
        point4.X() = std::move(tmp);
        root = square(std::move(root))%*gwstate.N;
        tmp = std::move(point1.X());
        point1.X() = tmp*root%*gwstate.N;
        file_proofpoint.children()[1]->write(point1);
        point1.X() = std::move(tmp);
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
        root = std::move(task.state()->X());

        point4.X() = std::move(point4.X())*root%*gwstate.N;
        file_proofpoint.children()[4]->write(point4);
        file_proofpoint.children()[2]->read(point2);
        point2.X() = std::move(point2.X())*root%*gwstate.N;
        file_proofpoint.children()[2]->write(point2);
        file_proofpoint.children()[1]->write(point1);
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
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());

    proof.points().clear();
    Fermat fermat2(Fermat::AUTO, input, params, logging, &proof);
    proof_build.points() = proof.points();

    fingerprint = input.fingerprint();
    file_cert = File("prst_cert", fingerprint);
    fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat2.a()) + "." + std::to_string(proof.points()[proof_count]));
    file_proofpoint = File("prst_proof", fingerprint);
    file_proofproduct = File("prst_prod", fingerprint);
    file_checkpoint = File("prst_c", fingerprint);
    file_recoverypoint = File("prst_r", fingerprint);

    proof.init_files(&file_proofpoint, &file_proofproduct, &file_cert);
    proof_build.init_files(&file_proofpoint, &file_proofproduct, &file_cert);

    try
    {
        fermat2.run(input, gwstate, file_checkpoint, file_recoverypoint, logging, &proof);

        FastExp task(power(input.gb()/3, 3)*power(input.gb(), input.n() - 3));
        task.init(&input, &gwstate, nullptr, &logging, fermat2.a());
        task.run();
        root = std::move(task.state()->X());

        file_proofpoint.children()[4]->read(point4);
        point4.X() = std::move(point4.X())*root%*gwstate.N;
        file_proofpoint.children()[4]->write(point4);
        root = square(std::move(root))%*gwstate.N;
        root = square(std::move(root))%*gwstate.N;
        file_proofpoint.children()[3]->read(point3);
        point3.X() = std::move(point3.X())*root%*gwstate.N;
        file_proofpoint.children()[3]->write(point3);
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
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());

    params.ProofPointsPerCheck = 1;
    proof.points().clear();
    Fermat fermat3(Fermat::AUTO, input, params, logging, &proof);
    proof_build.points() = proof.points();

    fingerprint = input.fingerprint();
    file_cert = File("prst_cert", fingerprint);
    fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat3.a()) + "." + std::to_string(proof.points()[proof_count]));
    file_proofpoint = File("prst_proof", fingerprint);
    file_proofproduct = File("prst_prod", fingerprint);
    file_checkpoint = File("prst_c", fingerprint);
    file_recoverypoint = File("prst_r", fingerprint);

    proof.init_files(&file_proofpoint, &file_proofproduct, &file_cert);
    proof_build.init_files(&file_proofpoint, &file_proofproduct, &file_cert);

    try
    {
        fermat3.run(input, gwstate, file_checkpoint, file_recoverypoint, logging, &proof);

        FastExp task((input.value() - 1)/3);
        task.init(&input, &gwstate, nullptr, &logging, 2);
        task.run();
        root = std::move(task.state()->X());

        file_proofpoint.children()[4]->read(point4);
        point4.X() = std::move(point4.X())*root%*gwstate.N;
        file_proofpoint.children()[4]->write(point4);
        file_proofpoint.children()[3]->read(point3);
        point3.X() = std::move(point3.X())*root%*gwstate.N;
        file_proofpoint.children()[3]->write(point3);
        root = square(std::move(root))%*gwstate.N;
        file_proofpoint.children()[2]->read(point2);
        point2.X() = std::move(point2.X())*root%*gwstate.N;
        file_proofpoint.children()[2]->write(point2);
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
