#ifndef _PROBABILISTIC_FILTER_H_
#define _PROBABILISTIC_FILTER_H_

#include <vector>
#include "../model/model.h"

// computes "optimal" parameters update based on measurements
class Probabilistic_Filter
{

protected:
    // forward function computing
    Model *m;

    // noise in measurements
    std::vector<double> R;

    // noise in parameters
    std::vector<double> S;

public:

    // get updates for parameter vector and covariance matrix
    virtual void get_update(double *mes, double *upd_vec, double *upd_cov) = 0;

    Probabilistic_Filter(Model *model) 
    : m(model), R(m->forward_size*m->forward_size, 0), S(m->num_pars*m->num_pars, 0) {};

    virtual ~Probabilistic_Filter() = default;
};
#endif