#ifndef _KALMAN_UNSCENTED_H_
#define _KALMAN_UNSCENTED_H_

// unscented kalman filter implementation

#include "../types/model/model.h"
#include "../types/filter/probabilistic_filter.h"
#include <cstring>

class kalman_unscented : public probabilistic_filter
{

private:

    // cross-covariances
    double P_y[m->forward_size*m->forward_size]
    double P_xy[m->num_pars*m->forward_size];             


    // matrix for computation purposes
    double P_aux[m->forward_size*m->forward_size];
    double P[m->num_pars*m->num_pars];
    double S_aux_1[m->num_pars*m->forward_size];    
    double P_y_aux[m->forward_size*m->forward_size];

    // kalman gain
    double K[m->num_pars*m->forward_size];

    // sigma points
    double x_sigma[2*m->num_pars*m->num_pars], 
    double y_sigma[2*m->num_pars*m->forward_size];      

    // averaging result
    double y_avg[m->forward_size];      
    
    // sigma points and response computation
    void compute_sigma();
   
    // covariance matrix calculation for inversion
    void get_matrices();

};
#endif