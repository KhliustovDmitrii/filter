#ifndef _SIMPLE_PARSER_H_
#define _SIMPLE_PARSER_H_

#include "../../../types/reader/reader.h"
// reads input string and parses to label and numeric data

class Simple_Parser : public Reader
{

public:
    // control codes:
    // 0  - succesfully read and parsed line from file
    // 1  - line is to be skipped
    // -1 - error during parsing
    virtual int read(std::vector<double> &dest, std::string &label) override;

    Simple_Parser(std::string fname, char s, char d, int ls, 
                  std::vector<std::string> dr) :
    Reader(fname),
    separator(s), decimal(d), label_size(ls), drop(dr) {};

private:

    std::string line;

    // supports separator consisting of characters of single types only
    // - I don't want to implement Aho-Korasik algorithm here
    char separator;            // separator between consecutive data values in file
    char decimal;              // the decimal . or , used in number representation
    int label_size;            // size of label in the beginning of line

    std::vector<std::string> drop;    // we drop lines containing these characters
                                      // - usually comments or missing measurements

};
#endif