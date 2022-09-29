
#include <deque>
#include <tuple>
#include "gwnum.h"
#include "cpuid.h"
#include "arithmetic.h"
#include "exception.h"
#include "inputnum.h"
#include "file.h"
#include "logging.h"
#include "task.h"
#include "params.h"
#include "fermat.h"
#include "proof.h"
#include "pocklington.h"
#include "testing.h"

#include "test.data"

using namespace arithmetic;

int testing_main(int argc, char *argv[])
{
    int i;
    Task::PROGRESS_TIME = 60;
    int log_level = Logging::LEVEL_WARNING;
    Params params;
    uint64_t maxMem = 0;
    std::string subset;

    for (i = 1; i < argc; i++)
        if (argv[i][0] == '-' && argv[i][1])
        {
            switch (argv[i][1])
            {
            case 't':
                if (argv[i][2] && isdigit(argv[i][2]))
                    params.thread_count = atoi(argv[i] + 2);
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    params.thread_count = atoi(argv[i]);
                }
                else
                    break;
                if (params.thread_count == 0 || params.thread_count > 64)
                    params.thread_count = 1;
                continue;

            case 'd':
                if (!argv[i][2])
                    log_level = Logging::LEVEL_INFO;
                else
                    break;
                continue;
            }

            if (i < argc - 1 && strcmp(argv[i], "-M") == 0)
            {
                i++;
                maxMem = InputNum::parse_numeral(argv[i]);
            }
            else if (i < argc - 1 && strcmp(argv[i], "-spin") == 0)
            {
                i++;
                params.spin_threads = atoi(argv[i]);
            }
            else if (i < argc - 1 && strcmp(argv[i], "-test") == 0)
            {
                i++;
                subset = argv[i];
            }
            else if (i < argc - 1 && strcmp(argv[i], "-check") == 0)
            {
                while (true)
                    if (i < argc - 1 && strcmp(argv[i + 1], "near") == 0)
                    {
                        i++;
                        params.CheckNear = true;
                        params.Check = false;
                    }
                    else if (i < argc - 1 && strcmp(argv[i + 1], "always") == 0)
                    {
                        i++;
                        params.CheckNear = false;
                        params.Check = true;
                    }
                    else if (i < argc - 1 && strcmp(argv[i + 1], "never") == 0)
                    {
                        i++;
                        params.CheckNear = false;
                        params.Check = false;
                    }
                    else if (i < argc - 1 && strcmp(argv[i + 1], "Gerbicz") == 0)
                    {
                        i++;
                        params.CheckGerbicz = true;
                    }
                    else
                        break;
            }
            else if (strcmp(argv[i], "-time") == 0)
            {
                while (true)
                    if (i < argc - 2 && strcmp(argv[i + 1], "write") == 0)
                    {
                        i += 2;
                        Task::DISK_WRITE_TIME = atoi(argv[i]);
                    }
                    else if (i < argc - 2 && strcmp(argv[i + 1], "progress") == 0)
                    {
                        i += 2;
                        Task::PROGRESS_TIME = atoi(argv[i]);
                    }
                    else
                        break;
            }
            else if (i < argc - 1 && strcmp(argv[i], "-log") == 0)
            {
                i++;
                if (strcmp(argv[i], "debug") == 0)
                    log_level = Logging::LEVEL_DEBUG;
                if (strcmp(argv[i], "info") == 0)
                    log_level = Logging::LEVEL_INFO;
                if (strcmp(argv[i], "warning") == 0)
                    log_level = Logging::LEVEL_WARNING;
                if (strcmp(argv[i], "error") == 0)
                    log_level = Logging::LEVEL_ERROR;
            }
        }
    if (subset.empty())
    {
        printf("Usage: PRST -test <subset> <options>\n");
        printf("Options: [-t <threads>] [-spin <threads>] [-log {debug | info | warning | error}] [-time [write <sec>] [progress <sec>]]\n");
        printf("\t-check [{near | always| never}] [Gerbicz] \n");
        printf("Subsets:\n");
        printf("\tall = 321plus + 321minus + b5plus + b5minus + gfn13 + special + error + prime\n");
        printf("\tslow = gfn13more + 100186b5minus + 109208b5plus\n");
        printf("\trandom\n");
        return 0;
    }

    TestLogging logging(log_level == Logging::LEVEL_ERROR ? Logging::LEVEL_ERROR : Logging::LEVEL_INFO);
    
    std::list<std::tuple<std::string,SubLogging,std::deque<Test>>> tests;
    auto add = [&](const std::string& subset) -> std::deque<Test>& { return std::get<2>(tests.emplace_back(subset, SubLogging(logging, log_level), std::deque<Test>())); };

    if (subset == "all" || subset == "321plus")
    {
        auto& cont = add("321plus");
        for (NTest* nTest = Test321Plus; nTest->n != 0; nTest++)
            cont.emplace_back(3, 2, *nTest, 1);
    }

    if (subset == "all" || subset == "321minus")
    {
        auto& cont = add("321minus");
        for (NTest* nTest = Test321Minus; nTest->n != 0; nTest++)
            cont.emplace_back(3, 2, *nTest, -1);
    }

    if (subset == "all" || subset == "b5plus")
    {
        auto& cont = add("b5plus");
        for (NTest* nTest = TestBase5Plus; nTest->n != 0; nTest++)
            cont.emplace_back(2, 5, *nTest, 1);
    }

    if (subset == "all" || subset == "b5minus")
    {
        auto& cont = add("b5minus");
        for (NTest* nTest = TestBase5Minus; nTest->n != 0; nTest++)
            cont.emplace_back(2, 5, *nTest, -1);
    }

    if (subset == "all" || subset == "gfn13")
    {
        auto& cont = add("gfn13");
        for (BTest* bTest = TestGFN13; bTest->b != 0; bTest++)
            cont.emplace_back(1, *bTest, 8192, 1);
    }

    if (subset == "all" || subset == "special")
    {
        auto& cont = add("special");
        for (KBNCTest* kbncTest = TestSpecial; kbncTest->n != 0; kbncTest++)
            cont.emplace_back(*kbncTest);
    }

    if (subset == "all" || subset == "error")
    {
        auto& cont = add("error");
    }
    
    if (subset == "all" || subset == "prime")
    {
        auto& cont = add("prime");
        for (KBNCTest* kbncTest = TestPrime; kbncTest->n != 0; kbncTest++)
            cont.emplace_back(*kbncTest);
    }
    
    if (subset == "slow" || subset == "gfn13more")
    {
        auto& cont = add("gfn13more");
        for (BTest* bTest = TestGFN13More; bTest->b != 0; bTest++)
            cont.emplace_back(1, *bTest, 8192, 1);
    }

    if (subset == "slow" || subset == "109208b5plus")
    {
        auto& cont = add("109208b5plus");
        for (KBNTest* kbnTest = Test109208Base5Plus; kbnTest->n != 0; kbnTest++)
            cont.emplace_back(*kbnTest, 1);
    }

    if (subset == "slow" || subset == "100186b5minus")
    {
        auto& cont = add("100186b5minus");
        for (KBNTest* kbnTest = Test100186Base5Minus; kbnTest->n != 0; kbnTest++)
            cont.emplace_back(*kbnTest, -1);
    }

    for (auto& subsetTests : tests)
    {
        for (auto& test : std::get<2>(subsetTests))
            std::get<1>(subsetTests).progress().add_stage(test.cost());
        logging.progress().add_stage(std::get<1>(subsetTests).progress().cost_total());
    }

    try
    {
        for (auto& subsetTests : tests)
        {
            logging.progress().update(0, 0);
            logging.warning("Running %s tests.\n", std::get<0>(subsetTests).data());
            if (std::get<0>(subsetTests) == "error")
            {
                SubLogging subLogging(std::get<1>(subsetTests), log_level > Logging::LEVEL_INFO ? Logging::LEVEL_ERROR + 1 : log_level);
                RootsTest(subLogging, params);
            }
            else
                for (auto& test : std::get<2>(subsetTests))
                {
                    std::get<1>(subsetTests).progress().update(0, 0);
                    SubLogging subLogging(std::get<1>(subsetTests), log_level > Logging::LEVEL_INFO ? Logging::LEVEL_ERROR : log_level);
                    test.run(subLogging, params);
                    std::get<1>(subsetTests).progress().next_stage();
                }
            logging.progress().next_stage();
        }
        logging.warning("All tests completed successfully.\n");
    }
    catch (const TaskAbortException&)
    {
        if (Task::abort_flag())
            logging.error("Test aborted.\n");
        else
            logging.error("Test failed.\n");
    }

    return 0;
}

void Test::run(Logging& logging, Params& global)
{
    Params params;
    params.Check = global.Check;
    params.CheckNear = global.CheckNear;
    params.CheckGerbicz = global.CheckGerbicz;
    if (input.b() <= 5) // tail mode
        params.ProofPointsPerCheck = 1;
    params.ProofSecuritySeed = "12345";
    params.RootOfUnityCheck = false;
    int proof_count = 16;

    Proof proof(Proof::SAVE, proof_count, input, params, logging);
    Fermat fermat(Fermat::AUTO, input, params, logging, &proof);

    uint32_t fingerprint = input.fingerprint();
    File file_cert("prst_cert", fingerprint);
    fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat.a()) + "." + std::to_string(proof.points()[proof_count]));
    File file_proofpoint("prst_proof", fingerprint);
    File file_proofproduct("prst_prod", fingerprint);
    File file_checkpoint("prst_c", fingerprint);
    File file_recoverypoint("prst_r", fingerprint);

    GWState gwstate;
    gwstate.thread_count = global.thread_count;
    gwstate.spin_threads = global.spin_threads;
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

        Proof proof_build(Proof::BUILD, proof_count, input, params, logging);
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

        Proof proof_cert(Proof::CERT, Proof::read_cert_power(file_cert), input, params, logging);
        proof_cert.run(input, gwstate, file_cert, file_checkpoint, file_recoverypoint, logging);
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

void RootsTest(Logging& logging, Params& global)
{
    Params params;
    params.Check = global.Check;
    params.CheckNear = global.CheckNear;
    params.CheckGerbicz = global.CheckGerbicz;
    params.RootOfUnityCheck = false;
    InputNum input;
    int proof_count = 4;
    Giant tmp, root;
    BaseExp::State point1;
    BaseExp::State point2;
    BaseExp::State point3;
    BaseExp::State point4;

    input.parse("3*2^353+1");

    Proof proof(Proof::SAVE, proof_count, input, params, logging);
    Fermat fermat(Fermat::AUTO, input, params, logging, &proof);
    Proof proof_build(Proof::BUILD, proof_count, input, params, logging);
    proof_build.points() = proof.points();

    uint32_t fingerprint = input.fingerprint();
    File file_cert("prst_cert", fingerprint);
    fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat.a()) + "." + std::to_string(proof.points()[proof_count]));
    File file_proofpoint("prst_proof", fingerprint);
    File file_proofproduct("prst_prod", fingerprint);
    File file_checkpoint("prst_c", fingerprint);
    File file_recoverypoint("prst_r", fingerprint);

    proof.init_files(&file_proofpoint, &file_proofproduct, &file_cert);
    proof_build.init_files(&file_proofpoint, &file_proofproduct, &file_cert);

    GWState gwstate;
    gwstate.thread_count = global.thread_count;
    gwstate.spin_threads = global.spin_threads;
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

        Proof proof_cert(Proof::CERT, Proof::read_cert_power(file_cert), input, params, logging);
        proof_cert.run(input, gwstate, file_cert, file_checkpoint, file_recoverypoint, logging);
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

        Proof proof_cert(Proof::CERT, Proof::read_cert_power(file_cert), input, params, logging);
        proof_cert.run(input, gwstate, file_cert, file_checkpoint, file_recoverypoint, logging);
        if (proof_cert.res64() != proof_build.res64())
        {
            logging.error("Attacking certificate mismatch.\n");
            throw TaskAbortException();
        }

        Proof proof_root(Proof::ROOT, proof_count, input, params, logging);
        try
        {
            proof_root.run(input, gwstate, logging, &fermat.result());
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
