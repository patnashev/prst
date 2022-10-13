#define NET_PRST_VERSION "0.7.0"
#define _SILENCE_CXX17_ALLOCATOR_VOID_DEPRECATION_WARNING

#include <stdio.h>
#include <stdlib.h>
#include "gwnum.h"
#include "cpuid.h"

#include "net.h"
#include "md5.h"
#include "task.h"
#include "exception.h"
#include "fermat.h"
#include "pocklington.h"
#include "proof.h"
#include "params.h"

using namespace restc_cpp;
using namespace arithmetic;

Writer* NetLogging::LoggingNetFile::get_writer()
{
    std::lock_guard<std::mutex> lock(net().upload_mutex());
    net().upload_cancel(this);
    return new Writer(std::move(_buffer));
}

void NetLogging::LoggingNetFile::on_upload()
{
    _buffer.swap(net().buffer());
    _buffer.clear();
}

void NetLogging::report(const std::string& message, int level)
{
    Logging::report(message, level);
    if (_net_level > level || level == LEVEL_RESULT)
        return;
    _file.write_text(message);
}

void NetLogging::report_param(const std::string& name, int value)
{
    if (name == "fft_len")
        _net.task()->fft_len = value;
    if (name == "a")
        _net.task()->a = value;
    if (name == "L")
        _net.task()->L = value;
    if (name == "L2")
        _net.task()->L2 = value;
    if (name == "M")
        _net.task()->M = value;
}

void NetLogging::report_param(const std::string& name, const std::string& value)
{
    if (name == "fft_desc")
        _net.task()->fft_desc = value;
}

void NetLogging::progress_save()
{
    _net.task()->progress = progress().progress_total();
    _net.task()->time = progress().time_total();
    _net.task()->time_op = progress().time_op()*1000;
}

void NetFile::on_upload()
{
    _uploading = true;
}

void NetFile::on_uploaded()
{
    _uploading = false;
    if (_free_buffer)
    {
        _free_buffer = false;
        std::vector<char>().swap(_buffer);
    }
}

File* NetFile::add_child(const std::string& name, uint32_t fingerprint)
{
    _children.emplace_back(new NetFile(_net_ctx, _filename + "." + name, fingerprint));
    return _children.back().get();
}

Writer* NetFile::get_writer()
{
    std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
    _free_buffer = false;
    _net_ctx.upload_cancel(this);
    if (_uploading)
    {
        _buffer.swap(_net_ctx.buffer());
        _uploading = false;
    }
    Writer* writer = new Writer(std::move(_buffer));
    writer->buffer().clear();
    return writer;
}

void NetFile::read_buffer()
{
    if (_free_buffer)
    {
        std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
        _free_buffer = false;
    }
    if (!_buffer.empty())
        return;

    boost::optional<std::string> md5;
    std::string str_data;

    // Run our example in a lambda co-routine
    auto done = _net_ctx.client()->ProcessWithPromiseT<bool>([&](Context& ctx) {
        // This is the co-routine, running in a worker-thread

        while (true)
            try
            {
                // Construct a request to the server
                auto reply = RequestBuilder(ctx)
                    .Get(_net_ctx.url() + "llr/" + _net_ctx.task_id() + "/" + filename())
                    .Argument("workerID", _net_ctx.worker_id())

                    // Send the request
                    .Execute();

                md5 = reply->GetHeader("MD5");
                str_data = reply->GetBodyAsString();

                return true;
            }
            catch (const HttpNotFoundException&) {
                //clog << "No file." << endl;
                str_data.clear();
                return true;
            }
            catch (const HttpForbiddenException&) {
                //clog << "No task." << endl;
                Task::abort();
                return false;
            }
            catch (const std::exception& ex) {
                std::clog << "File " << filename() << " download failed: " << ex.what() << std::endl;
                ctx.Sleep(boost::posix_time::microseconds(15000000));
                continue;
            }

        return true;
        });

    if (!done.get())
        return;

    if (hash && !str_data.empty() && md5)
    {
        char md5hash[33];
        md5_raw_input(md5hash, (unsigned char*)str_data.data(), (int)str_data.size());
        _md5hash = md5hash;
        if (md5.get() != _md5hash)
        {
            clear();
            return;
        }
    }

    _buffer.insert(_buffer.end(), str_data.begin(), str_data.end());
}

void NetFile::commit_writer(Writer& writer)
{
    std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
    _free_buffer = false;
    if (_uploading)
    {
        _buffer.swap(_net_ctx.buffer());
        _uploading = false;
    }
    if (hash)
        _md5hash = writer.hash_str();
    _buffer = std::move(writer.buffer());
    _net_ctx.upload(this);
}

void NetFile::free_buffer()
{
    std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
    if (_net_ctx.upload_queued(this))
    {
        _free_buffer = true;
        return;
    }
    if (_uploading)
    {
        _buffer.swap(_net_ctx.buffer());
        _uploading = false;
    }
    _free_buffer = false;
    std::vector<char>().swap(_buffer);
}

void NetFile::clear(bool recursive)
{
    std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
    _net_ctx.upload_cancel(this);
    if (_uploading)
    {
        _buffer.swap(_net_ctx.buffer());
        _uploading = false;
    }
    std::vector<char>().swap(_buffer);
    _md5hash.clear();
    _net_ctx.upload(this);
}

File* LLR2NetFile::add_child(const std::string& name, uint32_t fingerprint)
{
    _children.emplace_back(new LLR2NetFile(_net_ctx, _filename + "." + name, fingerprint, _type));
    _children.back()->hash = hash;
    return _children.back().get();
}

void LLR2NetFile::read_buffer()
{
    NetFile::read_buffer();

    if (_buffer.size() > 16 && *(uint32_t*)_buffer.data() == MAGIC_NUM && _buffer[4] == 2)
    {
        _buffer[4] = 4;
        _buffer[6] = _type;
        if (_type == BaseExp::State::TYPE)
            (*(uint32_t*)(_buffer.data() + 12))--;
    }
}

void LLR2NetFile::commit_writer(Writer& writer)
{
    if (writer.buffer().size() > 16)
    {
        writer.buffer()[4] = 2;
        writer.buffer()[6] = 0;
        if (_type == BaseExp::State::TYPE)
            (*(uint32_t*)(writer.buffer().data() + 12))++;
        writer.write((uint32_t)0);
        uint32_t checksum = 0;
        char* end = writer.buffer().data() + writer.buffer().size();
        for (char* data = writer.buffer().data() + 8; data < end; data += 4)
            checksum += *(uint32_t*)data;
        writer.write(checksum);
        for (int i = 0; i < 20; i++)
            writer.write((uint32_t)0);
    }

    NetFile::commit_writer(writer);
}

bool NetContext::upload_queued(NetFile* file)
{
    for (auto it = _upload_queue.begin(); it != _upload_queue.end(); it++)
        if (*it == file)
            return true;
    return false;
}

void NetContext::upload_cancel(NetFile* file)
{
    auto it = _upload_queue.begin();
    while (it != _upload_queue.end())
        if (*it == file)
            it = _upload_queue.erase(it);
        else
            it++;
}

void NetContext::upload_wait()
{
    std::future<void> localF;
    {
        std::lock_guard<std::mutex> lock(_upload_mutex);
        localF = std::move(_uploadF);
    }
    if (localF.valid())
        localF.get();
}

void NetContext::upload(NetFile* file)
{
    _upload_queue.push_back(file);
    if (_uploadF.valid() || _task->aborted)
        return;

    _uploadF = _putter->ProcessWithPromise([this](Context& ctx) {

        NetFile* file = nullptr;
        std::string put_url;
        std::string md5;
        char* data;
        size_t size;
        while (true)
        {
            if (file == nullptr)
            {
                std::lock_guard<std::mutex> lock(_upload_mutex);
                if (_upload_queue.empty() || _task->aborted)
                {
                    _uploadF = std::future<void>();
                    return;
                }
                file = _upload_queue.front();
                _upload_queue.pop_front();
                put_url = url() + "llr/" + task_id() + "/" + file->filename();
                data = file->buffer().data();
                size = file->buffer().size();
                md5 = file->md5hash();
                file->on_upload();
            }

            try
            {
                RequestBuilder(ctx)
                    .Put(put_url)
                    .Argument("md5", md5)
                    .Argument("workerID", worker_id())
                    .Argument("uptime", uptime())
                    .Argument("fft_desc", _task->fft_desc)
                    .Argument("fft_len", _task->fft_len)
                    .Argument("a", _task->a)
                    .Argument("L", _task->L)
                    .Argument("L2", _task->L2)
                    .Argument("M", _task->M)
                    .Argument("progress", std::to_string(_task->progress))
                    .Argument("time", std::to_string(_task->time))
                    .Argument("time_op", std::to_string(_task->time_op))
                    .Body(std::unique_ptr<RequestBody>(new RequestBodyData(data, size)))
                    .Execute();

                std::lock_guard<std::mutex> lock(_upload_mutex);
                file->on_uploaded();
                file = nullptr;
            }
            catch (const HttpAuthenticationException&) {
                std::clog << "Task timed out." << std::endl;
                _task->aborted = true;
                Task::abort();
                std::lock_guard<std::mutex> lock(_upload_mutex);
                file->on_uploaded();
                file = nullptr;
            }
            catch (const HttpForbiddenException&) {
                std::clog << "Task not found." << std::endl;
                _task->aborted = true;
                Task::abort();
                std::lock_guard<std::mutex> lock(_upload_mutex);
                file->on_uploaded();
                file = nullptr;
            }
            catch (const std::exception& ex) {
                std::clog << "Upload to " << put_url << " failed: " << ex.what() << std::endl;
                ctx.Sleep(boost::posix_time::microseconds(15000000));
            }
        }
    });
}

void NetContext::done()
{
    _putter->CloseWhenReady(true);
}


int net_main(int argc, char *argv[])
{
    int i;
    GWState gwstate;
    std::string url;
    std::string worker_id;
    int log_level = Logging::LEVEL_INFO;
    int net_log_level = Logging::LEVEL_WARNING;
    uint64_t maxMem = 2048*1048576ULL;
    int disk_write_time = Task::DISK_WRITE_TIME;

    for (i = 1; i < argc; i++)
        if (argv[i][0] == '-' && argv[i][1])
        {
            switch (argv[i][1])
            {
            case 't':
                if (argv[i][2] && isdigit(argv[i][2]))
                    gwstate.thread_count = atoi(argv[i] + 2);
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    gwstate.thread_count = atoi(argv[i]);
                }
                else
                    break;
                if (gwstate.thread_count == 0 || gwstate.thread_count > 64)
                    gwstate.thread_count = 1;
                continue;

            case 'M':
                if (argv[i][2] && isdigit(argv[i][2]))
                    maxMem = InputNum::parse_numeral(argv[i] + 2);
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    maxMem = InputNum::parse_numeral(argv[i]);
                }
                else
                    break;
                continue;

            case 'i':
                if (argv[i][2])
                    worker_id = argv[i] + 2;
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    worker_id = argv[i];
                }
                else
                    break;
                continue;
            }

            if (strcmp(argv[i], "-time") == 0)
            {
                while (true)
                    if (i < argc - 2 && strcmp(argv[i + 1], "write") == 0)
                    {
                        i += 2;
                        disk_write_time = atoi(argv[i]);
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
            else if (strcmp(argv[i], "-v") == 0)
            {
                printf("Net-PRST version " NET_PRST_VERSION ", Gwnum library version " GWNUM_VERSION "\n");
                return 0;
            }
        }
        else
        {
            url = argv[i];
        }
    if (url.empty() || url.find("http://") != 0 || worker_id.empty())
    {
        printf("Usage: PRST -net -i <WorkerID> http://<host>:<port>/api/\n");
        return 0;
    }

    NetContext net(url, worker_id, log_level, net_log_level);
    Logging& logging = net.logging();

	// Set the log-level to a reasonable value
	boost::log::core::get()->set_filter
	(
#ifdef _DEBUG
		boost::log::trivial::severity >= boost::log::trivial::warning
#else
		boost::log::trivial::severity >= boost::log::trivial::error
#endif // DEBUG
	);

	while (true)
	{
        if (Task::abort_flag())
            return 1;
        logging.set_prefix("");
        net.task().reset(new PRSTTask());

		// Run our example in a lambda co-routine
		auto done = net.client()->ProcessWithPromiseT<bool>([&](Context& ctx) {
			// This is the co-routine, running in a worker-thread

			try
			{

				// Construct a request to the server
				SerializeFromJson(*net.task(), RequestBuilder(ctx)
					.Post(net.url() + "llr/new")
					.Argument("workerID", net.worker_id())
					.Argument("uptime", net.uptime())
					.Argument("version", NET_PRST_VERSION)

					// Send the request
					.Execute()
				);
			}
			catch (const std::exception& ex) {
                std::clog << "Task acquisition failed: " << ex.what() << std::endl;
				return false;
			}

			return true;
		});
		if (!done.get())
		{
			std::this_thread::sleep_for(std::chrono::minutes(1));
			continue;
		}
        logging.info("%s\n", net.task_id().data());

        NetFile file_number(net, "number", 0);
        InputNum input;
        if (net.task()->n > 0)
            input.init(net.task()->sk, net.task()->sb, net.task()->n, net.task()->c);
        else if (!input.read(file_number))
        {
            logging.error("Number file is missing or corrupted.\n");
            std::this_thread::sleep_for(std::chrono::minutes(1));
            continue;
        }
        if (net.task()->options.find("FFT_Increment") != net.task()->options.end())
            gwstate.next_fft_count = std::stoi(net.task()->options["FFT_Increment"]);
        net.task()->a = net.task()->L = net.task()->L2 = net.task()->M = 0;
        int maxSize = (int)(maxMem/(gwnum_size(gwstate.gwdata())));

        logging.progress() = Progress();
        logging.progress().time_init(net.task()->time);
        if (net.task()->options.find("write_time") != net.task()->options.end())
            Task::DISK_WRITE_TIME = std::stoi(net.task()->options["write_time"]);
        else
            Task::DISK_WRITE_TIME = disk_write_time;
        
        Params params;
        bool supportLLR2 = true;
        if (net.task()->options.find("support") != net.task()->options.end())
            supportLLR2 = net.task()->options["support"] == "LLR2";
        int proof_op = Proof::NO_OP;
        if (net.task()->mode == "SavePoints")
            proof_op = Proof::SAVE;
        if (net.task()->mode == "BuildCert")
            proof_op = Proof::BUILD;
        if (net.task()->mode == "VerifyCert")
            proof_op = Proof::CERT;
        int proof_count = 16;
        if (net.task()->options.find("Gerbicz") != net.task()->options.end())
            params.CheckStrong = net.task()->options["Gerbicz"] == "1";
        if (net.task()->options.find("ProofCount") != net.task()->options.end())
            proof_count = std::stoi(net.task()->options["ProofCount"]);
        if (net.task()->options.find("PointsPerL2") != net.task()->options.end())
            params.StrongCount = proof_count/std::stoi(net.task()->options["PointsPerL2"]);
        if (net.task()->options.find("ProofName") != net.task()->options.end())
            params.ProofPointFilename = net.task()->options["ProofName"];
        if (net.task()->options.find("ProductName") != net.task()->options.end())
            params.ProofProductFilename = net.task()->options["ProductName"];

        std::list<std::unique_ptr<NetFile>> files;
        auto newFile = [&](const std::string& filename, uint32_t fingerprint, char type = BaseExp::State::TYPE)
        {
            if (supportLLR2)
                return files.emplace_back(new LLR2NetFile(net, filename, gwstate.fingerprint, type)).get();
            else
                return files.emplace_back(new NetFile(net, filename, fingerprint)).get();
        };
        uint32_t fingerprint = input.fingerprint();
        gwstate.fingerprint = fingerprint;
        File* file_cert = newFile(!params.ProofPointFilename.empty() ? params.ProofPointFilename + ".crt" : "proof.crt", fingerprint, Proof::Certificate::TYPE);
        std::unique_ptr<Proof> proof;
        if (proof_op != Proof::NO_OP)
            proof.reset(new Proof(proof_op, proof_count, input, params, *file_cert, logging));
        if (proof && net.task()->options.find("CachePoints") != net.task()->options.end())
            proof->set_cache_points(true);

        std::unique_ptr<Fermat> fermat;

        if (proof_op == Proof::CERT)
        {
        }
        else if (net.task()->type == "Pocklington")
        {
            fermat.reset(new Pocklington(input, params, logging, proof.get()));
        }
        else
            fermat.reset(new Fermat(Fermat::AUTO, input, params, logging, proof.get()));

        gwstate.maxmulbyconst = params.maxmulbyconst;
        input.setup(gwstate);
        logging.info("Using %s.\n", gwstate.fft_description.data());
        net.task()->fft_desc = gwstate.fft_description;
        net.task()->fft_len = gwstate.fft_length;

        try
        {
            if (proof_op == Proof::CERT)
            {
                fingerprint = File::unique_fingerprint(fingerprint, file_cert->filename());
                File* file_checkpoint = newFile("checkpoint", fingerprint);
                File* file_recoverypoint = newFile("recoverypoint", fingerprint);
                proof->run(input, gwstate, *file_checkpoint, *file_recoverypoint, logging);
            }
            else if (proof)
            {
                fingerprint = File::unique_fingerprint(fingerprint, std::to_string(fermat->a()) + "." + std::to_string(proof->points()[proof_count]));
                File* file_proofpoint = newFile(!params.ProofPointFilename.empty() ? params.ProofPointFilename : "proof", fingerprint);
                File* file_proofproduct = newFile(!params.ProofProductFilename.empty() ? params.ProofProductFilename : "prod", fingerprint, Proof::Product::TYPE);
                proof->init_files(file_proofpoint, file_proofproduct, file_cert);

                File* file_checkpoint = newFile("checkpoint", fingerprint);
                File* file_recoverypoint = newFile("recoverypoint", fingerprint);
                fermat->run(input, gwstate, *file_checkpoint, *file_recoverypoint, logging, proof.get());
            }
            else if (fermat)
            {
                File* file_checkpoint = newFile("checkpoint", fingerprint);
                File* file_recoverypoint = newFile("recoverypoint", fingerprint);
                fermat->run(input, gwstate, *file_checkpoint, *file_recoverypoint, logging, nullptr);
            }
        }
        catch (const TaskAbortException&)
        {
        }

        gwstate.done();

        net.upload_wait();
        if (net.task()->aborted)
        {
            Task::abort_reset();
            continue;
        }
        if (Task::abort_flag())
            return 1;

		// Run our example in a lambda co-routine
		auto doneRet = net.client()->ProcessWithPromise([&](Context& ctx) {
			// This is the co-routine, running in a worker-thread

			int i;
			for (i = 0; i < 10; i++)
				try
				{
					// Construct a request to the server
					RequestBuilder(ctx)
						.Post(net.url() + "llr/res/" + net.task_id())
						.Argument("workerID", net.worker_id())
                        .Argument("res", fermat->success() ? "prime" : fermat->res64())
                        .Argument("cert", proof ? proof->res64() : "")
                        .Argument("time", std::to_string(logging.progress().time_total()))
                        .Argument("version", NET_PRST_VERSION)

						// Send the request
						.Execute();
					break;
				}
				catch (const std::exception& ex) {
					std::clog << "Can't upload result: " << ex.what() << std::endl;
					ctx.Sleep(boost::posix_time::microseconds(i*5000000));
				}
		});

		doneRet.get();
	}

    net.done();

    return 0;
}

