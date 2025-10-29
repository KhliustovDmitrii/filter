#ifndef _KALMAN_EXTENDED_H_
#define _KALMAN_EXTENDED_H_

// extended kalman filter implementation

#include "../types/model/model.h"
#include "../types/filter/probabilistic_filter.h"
#include <cstring>

class kalman_extended : public probabilistic_filter
{

private:

    // jacobian matrix
    double Jacobian[m->num_pars*m->forward_size];

    // compute jacobian matrix
    void get_jacobian();

    // inner fuction for calculating update
    void proc(int i, double *upd_vec, double *upd_cov, double delta);

};
#endif