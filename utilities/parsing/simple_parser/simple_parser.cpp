#include "simple_parser.h"

int Simple_Parser::read(std::vector<double> &dest, std::string &label)
{
    if(!std::getline(file, line)) return -1;

    // drop conditions
    for(std::string str : drop)
        if(line.find(str) != std::string::npos) return 1;

    // cut off the label
    label = line.substr(0, label_size);

    // finite automata for string parsing
    int i, m_num, is_digit;
    std::string tmp;

    is_digit = 0;
    m_num = 0;
    for(i = label_size; i<line.size(); i++)
    {
        // boundary between numeric data pieces
        if(line[i] == separator)
        {
            if(is_digit == 1) // just finished reading number
            {
                dest[m_num] = std::stod(tmp);
                m_num++;
                tmp.clear();
            }
            is_digit = 0;
            continue;
        }

        // numeric character or decimal
        if(line[i] == decimal) tmp.push_back('.');
        else tmp.push_back(line[i]);

        is_digit = 1;
    }

    return 0;
}