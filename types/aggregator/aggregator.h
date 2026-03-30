#ifndef AGGREGATOR_H
#define AGGREGATOR_H

// class for aggregating several models
#include <vector>

#include "../model/model.h"

namespace filter::components
{
class Aggregator
{
public:
    virtual void aggregate(std::vector<double> weights) = 0;

    virtual ~Aggregator() = default;
};
}; // filter::components
#endif