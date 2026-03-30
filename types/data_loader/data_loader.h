#ifndef DATA_LOADER_H
#define DATA_LOADER_H

#include <vector>
#include "../model/model.h"

// class for splitting input data to model parameters and mesurements

namespace filter::io
{
class Data_Loader
{
public:

    // split data to model pars and measurements part
    virtual void load(const std::vector<double> & data) = 0;

    Data_Loader(int msize) :
    measurements(msize, 0.0) {};

    virtual ~Data_Loader() = default;

    void get_measurements(std::vector<double> &mes);

protected:
    std::vector<double> measurements;
};
}; // filter::io
#endif