#ifndef EQUATOR_AGGREGATOR_H
#define EQUATOR_AGGREGATOR_H

#include <vector>

#include "types/model/model.h"
#include "types/aggregator/aggregator.h"

#include "../EQUATOR_C/EQUATOR_C.h"

namespace filter::examples
{
// implementation of aggregating class for EQUATOR model
class EQUATOR_aggregator : public filter::components::Aggregator
{
public:
    virtual void aggregate(std::vector<double> weights) override;

    EQUATOR_aggregator(EQUATOR_C &m_out, std::vector<EQUATOR_C *> m_source) :
    out_model(m_out), source_models(m_source){};

private:
    EQUATOR_C &out_model;
    std::vector<EQUATOR_C *> source_models;
};
}; // filter::examples
#endif