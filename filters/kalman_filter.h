#ifndef _KALMAN_H_
#define _KALMAN_H_

#include "../models/model.h"
#include <cstring>

class kalman_filter
{
public:
   model<double, double> *m;
   
   double *R;
   double *S;
   double *Jacobian;
   double *lower_boundary;
   double *upper_boundary;
   
   void get_jacobian();
   void get_update_vector(double *upd_vec, double *measurements);
   void update_model(double *upd_vec, double *resp, double res);
   void proc(double *upd_vec, double delta, double *h, double r);
   
   kalman_filter(model<double ,double> *mod);
};



#endif
