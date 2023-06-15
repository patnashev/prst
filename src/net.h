#pragma once

#include <chrono>
#include <thread>
#include <deque>

#include "arithmetic.h"
#include "file.h"
#include "logging.h"
#include "inputnum.h"

#include <boost/lexical_cast.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>

#include "restc-cpp/restc-cpp.h"
#include "restc-cpp/SerializeJson.h"
#include "restc-cpp/RequestBuilder.h"


struct PRSTTask {
    std::string id;
    std::string type;
    std::string sk;
    std::string sb;
    int n = 0;
    int c = 0;
    int cyclotomic = 0;
    std::string proof;
    int count = 0;
    double time = 0.0;
    std::map<std::string, std::string> options;

    std::string fft_desc;
    int fft_len;
    int a, L, L2, M;
    double progress;
    double time_op;
    bool aborted = false;
};

BOOST_FUSION_ADAPT_STRUCT(
    PRSTTask,
    (auto, id)
    (auto, type)
    (auto, sk)
    (auto, sb)
    (auto, n)
    (auto, c)
    (auto, cyclotomic)
    (auto, proof)
    (auto, count)
    (auto, time)
    (auto, options)
)

class NetContext;

class NetFile : public File
{
public:
    NetFile(NetContext& net_ctx, const std::string& filename, uint32_t fingerprint) : File(filename, fingerprint), _net_ctx(net_ctx) { }

    File* add_child(const std::string& name, uint32_t fingerprint) override;

    using File::get_writer;
    Writer* get_writer() override;

    void read_buffer() override;
    void commit_writer(Writer& writer) override;
    void free_buffer() override;
    void clear(bool recursive = false) override;

    virtual void on_upload();
    virtual void on_uploaded();

    NetContext& net() { return _net_ctx; }
    std::string& md5hash() { return _md5hash; }

protected:
    NetContext& _net_ctx;
    std::string _md5hash;
    bool _uploading = false;
    bool _free_buffer = false;
};

class LLR2NetFile : public NetFile
{
public:
    LLR2NetFile(NetContext& net_ctx, const std::string& filename, uint32_t fingerprint, char type) : NetFile(net_ctx, filename, fingerprint), _type(type) { }

    File* add_child(const std::string& name, uint32_t fingerprint) override;
    void read_buffer() override;
    void commit_writer(Writer& writer) override;

protected:
    char _type;
};

class NetLogging : public Logging
{
public:
    class LoggingNetFile : public NetFile
    {
    public:
        LoggingNetFile(NetContext& net) : NetFile(net, "stderr", 0) { }

        Writer* get_writer() override;
        void on_upload() override;
    };

public:
    NetLogging(int level, int net_level, NetContext& net) : Logging(level), _net_level(net_level), _net(net), _file(net) { }

    virtual void report(const std::string& message, int level) override;
    virtual void report_param(const std::string& name, int value) override;
    virtual void report_param(const std::string& name, const std::string& value) override;
    virtual void report_progress() override;
    virtual void progress_save() override;

private:
    int _net_level;
    NetContext& _net;
    LoggingNetFile _file;
};

class NetContext
{
public:
    NetContext(std::string& url, std::string& worker_id, int log_level, int net_log_level) : _url(url), _worker_id(worker_id), _logging(log_level, net_log_level, *this), _start_time(std::chrono::system_clock::now())
    {
        _client = restc_cpp::RestClient::Create();

        restc_cpp::Request::Properties properties;
        properties.headers["Content-Type"] = "application/octet-stream";
        _putter = restc_cpp::RestClient::Create(properties);
    }

    void upload(NetFile* file);
    bool upload_queued(NetFile* file);
    void upload_cancel(NetFile* file);
    void upload_wait();
    void done();

    std::string& url() { return _url; }
    std::string& worker_id() { return _worker_id; }
    std::string& task_id() { return _task->id; }
    std::chrono::seconds::rep uptime() { return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - _start_time).count(); }

    NetLogging& logging() { return _logging; }
    restc_cpp::RestClient* client() { return _client.get(); }
    std::unique_ptr<PRSTTask>& task() { return _task; }

    std::mutex& upload_mutex() { return _upload_mutex; }
    std::vector<char>& buffer() { return _upload_buffer; }

private:
    void upload_written();

private:
    std::string _url;
    std::string _worker_id;
    NetLogging _logging;
    std::chrono::system_clock::time_point _start_time;

    std::unique_ptr<restc_cpp::RestClient> _client;
    std::unique_ptr<restc_cpp::RestClient> _putter;

    std::unique_ptr<PRSTTask> _task;

    std::future<void> _uploadF;
    std::deque<NetFile*> _upload_queue;
    std::mutex _upload_mutex;
    std::vector<char> _upload_buffer;
};

class RequestBodyData : public restc_cpp::RequestBody
{
public:
    RequestBodyData(const char *buffer, size_t count) : m_buffer{ buffer }, m_count{ count } { }
    Type GetType() const noexcept override { return Type::FIXED_SIZE; }
    std::uint64_t GetFixedSize() const override { return m_count; }
    void Reset() override { m_eof = false; }
    std::string GetCopyOfData() const override { return std::string(m_buffer, m_buffer + m_count); }

    bool GetData(restc_cpp::write_buffers_t & buffers) override
    {
        if (m_eof)
            return false;
        buffers.push_back({ m_buffer, m_count });
        m_eof = true;
        return true;
    }

private:
    const char *m_buffer;
    size_t m_count;
    bool m_eof = false;
};
