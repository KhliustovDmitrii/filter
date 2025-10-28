#ifndef _KALMAN_UNSCENTED_H_
#define _KALMAN_UNSCENTED_H_

#include "../models/model.h"
#include <cstring>

class kalman_unscented
{
public:
   model<double, double> *m;
   
   double *R;                      // measurement error
   double *S;                      // square root of parameter error
   double *P_y, *P_xy;             // cross-covariances
   double *P_aux, *P, *S_aux_1;    // matrix for computation purposes
   double *P_y_aux;
   double *K;                      // kalman gain
   double *x_sigma, *y_sigma;      // sigma points
   double *y_avg;                  // averaging result
   double *lower_boundary;         // parameter boundaries
   double *upper_boundary;
   
   // sigma points and response computation
   void compute_sigma();
   
   // covariance matrix calculation for inversion
   void get_matrices();
   
   // model update
   void get_update_vector(double *upd_vec, double *upd_cov, double *measurements);
   void update_cov(double *upd_cov, double factor);
   
   void invert_matrix(double *dest, double *source, int dim);
   void matrix_square_root(double *dest, double *source, int dim);
   
   kalman_unscented(model<double, double> *mod);
};



#endif
