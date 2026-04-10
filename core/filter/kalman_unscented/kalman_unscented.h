#ifndef KALMAN_UNSCENTED_H
#define KALMAN_UNSCENTED_H

// unscented kalman filter implementation

#include "types/model/model.h"
#include "types/filter/probabilistic_filter.h"
#include <cstring>

namespace filter
{
class Kalman_Unscented : public Probabilistic_Filter
{

private:

    struct UKF_Workspace : public Filter_Workspace
    {
        // cross-covariances
        std::vector<double> P_y, P_xy;

        // matrix for computation purposes
        std::vector<double> P_aux, P, S_aux_1, P_y_aux;

        // kalman gain
        std::vector<double> K;

        // sigma points
        std::vector<std::vector<double>> x_sigma, y_sigma;      

        // averaging result
        std::vector<double> y_avg;     
    };
    
    // sigma points and response computation
    void compute_sigma(Filter_Workspace &ws) const;
   
    // covariance matrix calculation for inversion
    void get_matrices(Filter_Workspace &ws) const;

public:

    void get_update(std::vector<double> &mes,
                    std::vector<double> &upd_vec, 
                    std::vector<double> &upd_cov,
                    Filter_Workspace &ws) const override;

    // update covariance matrix
    void update_covariance(Filter_Workspace &ws, std::vector<double> &upd_cov, double factor);

    Kalman_Unscented(Model &model) : Probabilistic_Filter(model) {};

    std::unique_ptr<Filter_Workspace> allocate_workspace() const override;
};
}; // filter
#endif