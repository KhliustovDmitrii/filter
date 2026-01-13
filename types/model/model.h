#ifndef _MODEL_H_
#define _MODEL_H_

#include <vector>

// computes response based on params
class Model
{

protected:

    // model parameters
    std::vector<double> params_;

public:

    // number of model parameters
    int num_pars;

    // dimensionality of response
    int forward_size;

    // forward function
    virtual void response(std::vector<double> &resp_arr) = 0;

    // proximity in response between model and measurements
    virtual double residual(std::vector<double> &mes, std::vector<double> &resp) = 0;

    // getters and setters
    virtual void set_param(int ind, double val) = 0;
    virtual double get_param(int ind) = 0;

    Model(int num_params = 0, int response_size = 0) 
        : num_pars(num_params), forward_size(response_size), params_(num_params, 0.0) {};

    virtual ~Model() = default;
};
#endif