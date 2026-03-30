#ifndef DECAY_UPDATER_H
#define DECAY_UPDATER_H
#include "types/updater/updater.h"
#include "types/model/model.h"

namespace filter
{
// updates model parameters with decaying step
class Decay_Updater : Updater
{
private:

    // number of decay steps
    int num_steps;

    // decay factor
    double factor;
public:
    Decay_Updater(Model &model, int ns, double fac) : 
    Updater(model),
    num_steps(ns),
    factor(fac){};

    double update(std::vector<double> &upd_vec, std::vector<double> &mes);
};
}; // filter
#endif