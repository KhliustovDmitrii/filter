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

    struct EKF_Workspace : public Filter_Workspace
    {
        std::vector<double> Jacobian;
        std::vector<double> e;
        std::vector<double> f;
    };

    // compute jacobian matrix
    void get_jacobian(EKF_Workspace &ws) const;

    // inner fuction for calculating update
    void proc(EKF_Workspace &ws, int k, 
              std::vector<double> &upd_vec, std::vector<double> &upd_cov, double delta) const;


public:
    Kalman_Extended(Model &model) 
    : Probabilistic_Filter(model){};

    std::unique_ptr<Filter_Workspace> allocate_workspace() const override;

    void get_update(std::vector<double> &mes,
                    std::vector<double> &upd_vec, 
                    std::vector<double> &upd_cov,
                    Filter_Workspace &ws) const override;
};
}; // filter
#endif