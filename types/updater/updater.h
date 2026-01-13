#ifndef _UPDATER_H_
#define _UPDATER_H_

#include "../model/model.h"

// updates model parameters
class Updater
{
    
protected:
    Model *m;

public:

    virtual double update(std::vector<double> &upd_vec, std::vector<double> &mes) = 0;

    Updater(Model *model) : m(model){};
    virtual ~Updater() = default;
};
#endif