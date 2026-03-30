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

    // cross-covariances
    mutable std::vector<double> P_y, P_xy;

    // matrix for computation purposes
    mutable std::vector<double> P_aux, P, S_aux_1, P_y_aux;

    // kalman gain
    mutable std::vector<double> K;

    // sigma points
    mutable std::vector<std::vector<double>> x_sigma, y_sigma;      

    // averaging result
    mutable std::vector<double> y_avg;     
    
    // sigma points and response computation
    void compute_sigma() const;
   
    // covariance matrix calculation for inversion
    void get_matrices() const;

public:

    void get_update(std::vector<double> &mes,
                    std::vector<double> &upd_vec, 
                    std::vector<double> &upd_cov) const override;

    // update covariance matrix
    void update_covariance(std::vector<double> &upd_cov, double factor);

    Kalman_Unscented(Model &model) : Probabilistic_Filter(model),
    P_y(model.forward_size*model.forward_size, 0),
    P_xy(m.num_pars*m.forward_size, 0),
    P_aux(m.forward_size*m.forward_size, 0),
    P(m.num_pars*m.num_pars, 0),
    S_aux_1(m.num_pars*m.forward_size, 0),
    P_y_aux(m.forward_size*m.forward_size),
    K(m.num_pars*m.forward_size, 0),
    x_sigma(2*m.num_pars, std::vector<double>(m.num_pars, 0)),
    y_sigma(2*m.num_pars, std::vector<double>(m.forward_size, 0)),
    y_avg(m.forward_size, 0) {};
};
}; // filter
#endif