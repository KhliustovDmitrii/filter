#ifndef _UNSCENTED_H_
#define _UNSCENTED_H_

#include "model.h"
#include <cstring>

/* class implementing unscented kalman filter algorithm */
class kalman_unscented
{
public:
   model *m;
   
   double *R;                      // measurement error
   double *S;                      // square root of parameter error
   
   // model update
   void get_update_vector(double *upd_vec, double *upd_cov, double *measurements);
   
   kalman_unscented(model *mod);
   
private:
   double *P_y, *P_xy;             // cross-covariances
   double *P_y_inv, *P, *S_aux;    // matrix for computation purposes
   double *P_y_aux;
   double *K;                      // kalman gain
   double *x_sigma, *y_sigma;      // sigma points
   double *y_avg;                  // averaging result
   
   // sigma points and response computation
   void compute_sigma();
   
   // covariance matrix calculation for inversion
   void get_matrices();
};



#endif
