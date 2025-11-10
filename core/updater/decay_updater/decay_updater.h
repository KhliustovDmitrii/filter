#ifndef _DECAY_UPDATER_H_
#define _DECAY_UPDATER_H_
#include "../types/updater/updater.h"
#include "../types/model/model.h"
// updates model parameters with decaying step
class Decay_Updater : Updater
{
private:

    // number of decay steps
    int num_steps;

    // decay factor
    double factor;
public:
    Decay_Updater(Model *model, int ns, int fac) : 
    Updater(model),
    num_steps(ns),
    factor(fac){};

    double update(double *upd_vec, double *mes);
};
#endif