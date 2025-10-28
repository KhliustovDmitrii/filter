#ifndef _COIL_MODEL_H_
#define _COIL_MODEL_H_

#include <cmath>
#include <iostream>
#include "model.h"

class coil_model : public model<double, double>
{
public:

   double coil_pos_1[3]; // first coil position
   double coil_pos_2[3]; // second coil position
   
   double mom_1[3];      // first coil moment
   double mom_2[3];      // second coil moment
   
   // auxillary variables for computations speedup
   double r_1[3];        // pendulum position w r t first coil
   double r_2[3];        // pendulum position w r t second coil
   double O[9];          // orientation matrix for receiver
   double tmp_vec[3];
   
   virtual void response(double *resp_arr);
   void quat_normalize(); // normalize quaternion to 1 norm
   
   coil_model()
   {
      int i;
      num_pars = 6;
      forward_size = 6;
      
      params = new double[num_pars];
      
      for(i=0; i<num_pars; i++) // pendulum positioned in the center of coordinate system
         params[i] = 0;
   }
};



#endif
