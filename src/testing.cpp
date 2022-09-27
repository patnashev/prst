
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
    int thread_count = 1;
    uint64_t maxMem = 0;
    std::string subset;

    for (i = 1; i < argc; i++)
        if (argv[i][0] == '-' && argv[i][1])
        {
            switch (argv[i][1])
            {
            case 't':
                if (argv[i][2] && isdigit(argv[i][2]))
                    thread_count = atoi(argv[i] + 2);
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    thread_count = atoi(argv[i]);
                }
                else
                    break;
                if (thread_count == 0 || thread_count > 64)
                    thread_count = 1;
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
            else if (i < argc - 1 && strcmp(argv[i], "-test") == 0)
            {
                i++;
                subset = argv[i];
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
        printf("Options: [-t <threads>] [-log {debug | info | warning | error}] [-time [write <sec>] [progress <sec>]]\n");
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
        for (KNTest* knTest = Test321Plus; knTest->n != 0; knTest++)
            cont.emplace_back(*knTest, 2, 1);
    }

    if (subset == "all" || subset == "321minus")
    {
        auto& cont = add("321minus");
        for (KNTest* knTest = Test321Minus; knTest->n != 0; knTest++)
            cont.emplace_back(*knTest, 2, -1);
    }

    if (subset == "all" || subset == "b5plus")
    {
        auto& cont = add("b5plus");
        for (KBNTest* kbnTest = TestBase5Plus; kbnTest->n != 0; kbnTest++)
            cont.emplace_back(*kbnTest, 1);
    }

    if (subset == "all" || subset == "b5minus")
    {
        auto& cont = add("b5minus");
        for (KBNTest* kbnTest = TestBase5Minus; kbnTest->n != 0; kbnTest++)
            cont.emplace_back(*kbnTest, -1);
    }

    if (subset == "all" || subset == "gfn13")
    {
        auto& cont = add("gfn13");
        for (KBNTest* kbnTest = TestGFN13; kbnTest->n != 0; kbnTest++)
            cont.emplace_back(*kbnTest, 1);
    }

    if (subset == "all" || subset == "special")
    {
        auto& cont = add("special");
        for (KBNCTest* kbncTest = TestSpecial; kbncTest->n != 0; kbncTest++)
            cont.emplace_back(*kbncTest);
    }

    if (subset == "all" || subset == "error")
    {
        //logging.error("Running tests with errors.\n");
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
        for (KBNTest* kbnTest = TestGFN13More; kbnTest->n != 0; kbnTest++)
            cont.emplace_back(*kbnTest, 1);
    }

    if (subset == "slow" || subset == "100186b5minus")
    {
        auto& cont = add("100186b5minus");
        for (KBNTest* kbnTest = Test100186Base5Minus; kbnTest->n != 0; kbnTest++)
            cont.emplace_back(*kbnTest, -1);
    }

    if (subset == "slow" || subset == "109208b5plus")
    {
        auto& cont = add("109208b5plus");
        for (KBNTest* kbnTest = Test109208Base5Plus; kbnTest->n != 0; kbnTest++)
            cont.emplace_back(*kbnTest, 1);
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
            logging.error("Running %s tests.\n", std::get<0>(subsetTests).data());
            for (auto& test : std::get<2>(subsetTests))
            {
                std::get<1>(subsetTests).progress().update(0, 0);
                SubLogging subLogging(std::get<1>(subsetTests), log_level > Logging::LEVEL_INFO ? Logging::LEVEL_ERROR : log_level);
                test.run(subLogging, thread_count);
                std::get<1>(subsetTests).progress().next_stage();
            }
        }
        logging.error("All tests completed successfully.\n");
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

void Test::run(Logging& logging, int thread_count)
{
    Params params;
    int proof_count = 16;
    params.ProofPointsPerCheck = 1;
    params.ProofSecuritySeed = "12345";

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
    gwstate.thread_count = thread_count;
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());

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
        file_cert.clear();
        file_proofpoint.clear(true);
        file_proofproduct.clear(true);
        file_checkpoint.clear(true);
        file_recoverypoint.clear(true);
        gwstate.done();
        throw;
    }

    gwstate.done();
}
