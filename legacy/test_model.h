#ifndef _TEST_MODEL_H_
#define _TEST_MODEL_H_


#include "model.h"

class test_model : public model<double, double>
{
public:
   void response(double *resp_arr)
   {
      resp_arr[0] = params[0]*3;
      resp_arr[1] = params[1];
      resp_arr[2] = params[0] + params[1] + 10;
   }
   
   test_model()
   {
      num_pars = 2;
      forward_size = 3;
      
      params = new double[num_pars];
      
      params[0] = 5;
      params[1] = -6;
   }
};

#endif
