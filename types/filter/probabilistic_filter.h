#ifndef _PROBABILISTIC_FILTER_H_
#define _PROBABILISTIC_FILTER_H_

#include "../model/model.h"

class probabilistic_filter
{
public:
    // forward function computing
    model *m;

    // noise in measurements
    double R[m->forward_size*m->forward_size];

    // noise in parameters
    double S[m->num_pars*m->num_pars];

    // get updates for parameter vector and covariance matrix
    virtual void get_update(double *mes, double *upd_vec, double *upd_cov) = 0;
}
#endif