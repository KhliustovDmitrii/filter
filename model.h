#ifndef _MODEL_H_
#define _MODEL_H_

/* general interface of model for solving forward problem */
class model
{
public:
   int num_pars;                              // number of parameters
   int forward_size;                          // dimensionality of output
   double *params;                            // model parameters
   virtual void response(double *resp_arr);   // forward problem solution
};

#endif
