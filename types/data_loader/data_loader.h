#ifndef _DATA_LOADER_H_
#define _DATA_LOADER_H_

#include <vector>
#include "../model/model.h"

// class for splitting input data to model parameters and mesurements

class Data_Loader
{
public:

    // split data to model pars and measurements part
    virtual void load(const std::vector<double> & data) = 0;

    Data_Loader(int msize) :
    measurements(msize, 0.0) {};

    virtual ~Data_Loader() = default;

protected:
    std::vector<double> measurements;
};

#endif