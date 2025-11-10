#ifndef _KALMAN_EXTENDED_H_
#define _KALMAN_EXTENDED_H_

#include <cstring>
#include <vector>
#include <cmath>
#include "../types/model/model.h"
#include "../types/filter/probabilistic_filter.h"

// extended Kalman filter implementation
class Kalman_Extended : public Probabilistic_Filter
{

private:

    // jacobian matrix
    std::vector<double> Jacobian;

    // compute jacobian matrix
    void get_jacobian();

    // inner fuction for calculating update
    void proc(int i, double *upd_vec, double *upd_cov, double delta);

public:
    Kalman_Extended(Model *model) 
    : Probabilistic_Filter(model), Jacobian(model->num_pars*model->forward_size, 0){};

    void get_update(double *mes, double *upd_vec, double *upd_cov);


};
#endif