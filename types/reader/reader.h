#ifndef _READER_H_
#define _READER_H_

#include <vector>
#include <fstream>

// reads data line by line, processes, stores into vector
class Reader
{
public:

    // interface method for reading next data line and storing parsed data in dest
    // returns control code describing outcome (e.g., success, no lines left etc.)
    // label saves non-numeric data needed to identify the line in file
    virtual int read(std::vector<double> &dest, std::string &label) = 0;

    Reader(std::string filename) : file(filename) 
    {
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
    };

    virtual ~Reader() = default;

protected:

    // source data file
    std::ifstream file;

};

#endif