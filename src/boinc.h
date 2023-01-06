#pragma once

int boinc_main(int argc, char *argv[]);

class BoincLogging : public Logging
{
public:
    BoincLogging() : Logging(LEVEL_INFO) { }

    void report(const std::string& message, int level) override;
    void report_progress() override;
    bool state_save_flag() override;
    void state_save() override;
    void progress_save() override;
    void heartbeat() override;
};
