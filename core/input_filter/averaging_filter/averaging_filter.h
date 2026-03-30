#ifndef AVERAGING_FILTER_H
#define AVERAGING_FILTER_H

#include <vector>
#include "types/input_filter/input_filter.h"

namespace filter::components
{
// simple input data filter
// averages specified number of input vectors
class Averaging_Filter : public Input_Filter
{
private:
    int window_size; // how many consecutive measurements to aggregate
    int agg_state; // how many measurements we have seen yet

    std::vector<double> mes_state; // sum and aggregation results

public:

    Averaging_Filter(int window_size_, int mes_size_) : 
    window_size(window_size_),
    mes_state(mes_size_, 0),
    agg_state(0){};

    // return codes:
    // -1   error adding data piece
    // 0    succesfully added data piece, aggregation result ready in output
    // 1    succesfully added, have more data to read (agg_state < window_size)
    virtual int add_data(std::vector<double> &data, std::vector<double> &output) override;
};
}; // filter::components
#endif