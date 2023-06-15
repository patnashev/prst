#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "gwnum.h"
#include "file.h"
#include "md5.h"
#include "inputnum.h"
#include "task.h"
#include "support.h"
#include "exp.h"
#include "proof.h"

File* LLR2File::add_child(const std::string& name, uint32_t fingerprint)
{
    _children.emplace_back(new LLR2File(_filename + "." + name, fingerprint, _type));
    _children.back()->hash = hash;
    return _children.back().get();
}

void LLR2File::read_buffer()
{
    File::read_buffer();

    if (_buffer.size() > 16 && *(uint32_t*)_buffer.data() == MAGIC_NUM && _buffer[4] == 2)
    {
        _buffer[4] = 4;
        _buffer[6] = _type;
        if (_type == BaseExp::StateValue::TYPE)
            (*(uint32_t*)(_buffer.data() + 12))--;
    }
}

void LLR2File::commit_writer(Writer& writer)
{
    if (writer.buffer().size() > 16)
    {
        writer.buffer()[4] = 2;
        writer.buffer()[6] = 0;
        if (_type == BaseExp::StateValue::TYPE)
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

    File::commit_writer(writer);
}
