#ifndef KALMAN_EXTENDED_H
#define KALMAN_EXTENDED_H

#include <cstring>
#include <vector>
#include <cmath>
#include "types/model/model.h"
#include "types/filter/probabilistic_filter.h"

namespace filter
{
// extended Kalman filter implementation
class Kalman_Extended : public Probabilistic_Filter
{

private:

    // jacobian matrix
    mutable std::vector<double> Jacobian;

    // compute jacobian matrix
    void get_jacobian() const;

    // inner fuction for calculating update
    void proc(int k, std::vector<double> &upd_vec, std::vector<double> &upd_cov, double delta) const;

public:
    Kalman_Extended(Model &model) 
    : Probabilistic_Filter(model), Jacobian(model.num_pars*model.forward_size, 0){};

    void get_update(std::vector<double> &mes,
                    std::vector<double> &upd_vec, 
                    std::vector<double> &upd_cov) const override;
};
}; // filter
#endif