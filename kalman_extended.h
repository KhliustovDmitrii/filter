#ifndef _EXTENDED_H_
#define _EXTENDED_H_

#include "model.h"
#include <cstring>

/* class implementing extended kalman filter algorithm */
class kalman_extended
{
public:
   model *m;                     // forward problem model
   
   double *R;                    // measurement noise covariance
   double *S;                    // parameter noise covariance
   
   // produce parameter update according to algorithm
   // TODO: addd posterior covariance supports
   void get_update_vector(double *upd_vec, double *measurements
                          int *valid_mes);
                          
   kalman_extended(model *mod);
   
private:
   double *Jacobian;             // Jacobi matrix
   double *f, *e;                // helper arrays
   
   // helper function for measurement processing
   void proc(double *upd_vec, double delta, double *h, double r);
};

#endif
