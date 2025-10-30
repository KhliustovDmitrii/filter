#ifndef _MODEL_H_
#define _MODEL_H_
class model
{
public:

    // number of model parameters
    int num_pars;

    // dimensionality of response
    int forward_size;

    // model parameters
    double params[num_pars];

    // forward function
    virtual void response(double *resp_arr) = 0;

    // update model parameters
    virtual void update(double *upd_vec) = 0;

    // proximity in response between model and measurements
    virtual double residual(double *mes, double *resp) = 0;

    virtual void set_param(int ind, double val) = 0;
    virtual double get_param(int ind) = 0;
}
#endif