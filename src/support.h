#pragma once

#include "file.h"

class LLR2File : public File
{
public:
    LLR2File(const std::string& filename, uint32_t fingerprint, char type) : File(filename, fingerprint), _type(type){ }

    File* add_child(const std::string& name, uint32_t fingerprint) override;
    void read_buffer() override;
    void commit_writer(Writer& writer) override;

protected:
    char _type;
};
