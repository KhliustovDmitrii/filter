#ifndef _INPUT_FILTER_H_
#define _INPUT_FILTER_H_

#include <vector>

// reads incoming data vector, processes, [outputs filtered data vector]
class Input_Filter
{
public:
    // get next data piece in data, [write filtering result to output], return control code
    // e.g. filtering result ready, more data needed, error processing
    virtual int add_data(std::vector<double> &data, std::vector<double> &output) = 0;

    virtual ~Input_Filter() = default;
};
#endif