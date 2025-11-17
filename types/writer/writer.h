#ifndef _WRITER_H_
#define _WRITER_H_

#include <vector>
#include <fstream>

#include "../model/model.h"

// stores model state to file
class Writer
{
public:

    virtual int write() = 0;

    Writer(std::string filename) : file(filename) 
    {
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
    };

    virtual ~Writer() = default;

protected:

    // file to output to
    std::ofstream file;

    Model *m;

};

#endif