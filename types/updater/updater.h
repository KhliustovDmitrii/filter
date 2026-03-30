#ifndef UPDATER_H
#define UPDATER_H

#include "../model/model.h"

namespace filter
{
// updates model parameters
class Updater
{
    
protected:
    Model &m;

public:

    virtual double update(std::vector<double> &upd_vec, std::vector<double> &mes) = 0;

    Updater(Model &model) : m(model){};
    virtual ~Updater() = default;
};
}; // filter
#endif