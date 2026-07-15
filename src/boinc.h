#pragma once

int boinc_main(int argc, char *argv[]);

class BoincLogging : public Logging
{
public:
    BoincLogging() : Logging(LEVEL_INFO) { _file_prime.clear(); _file_factor.clear(); }

    void report(const std::string& message, int level) override;
    void report_progress() override;
    bool state_save_flag() override;
    void progress_save() override;
    void heartbeat() override;

    bool send_trickle_messages = true;
};
