#ifndef _AGGREGATOR_H_
#define _AGGREGATOR_H_

// class for aggregating several models
#include <vector>

#include "../model/model.h"

class Aggregator
{
public:
    virtual void aggregate(std::vector<double> weights) = 0;

    virtual ~Aggregator() = default;
};

#endif