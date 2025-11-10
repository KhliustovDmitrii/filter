#ifndef _UPDATER_H_
#define _UPDATER_H_

#include "../model/model.h"

// updates model parameters
class Updater
{
    
protected:
    Model *m;

public:

    virtual double update(double *upd_vec, double *mes) = 0;

    Updater(Model *model) : m(model){};
    virtual ~Updater() = default;
};
#endif