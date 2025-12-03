#ifndef _POLYNOMIAL_MODEL_H_
#define _POLYNOMIAL_MODEL_H_

#include <vector>
#include "types/model/model.h"

// model computes values of polynomial in specified ponts
// residual is calculated in Lp norm
class Polynomial_Model : public Model
{
public:
    Polynomial_Model(std::vector<double> pts, std::vector<double> coeff, double norm) : 
    Model(), points(pts), norm_p(norm)
    {
        num_pars = coeff.size();
        forward_size = pts.size();

        params_ = coeff;
    };

    // forward function
    virtual void response(std::vector<double> &resp_arr) override;

    // proximity in response between model and measurements
    virtual double residual(std::vector<double> &mes, std::vector<double> &resp) override;

    // getters and setters
    virtual void set_param(int ind, double val) override;
    virtual double get_param(int ind) override;

private:
    std::vector<double> points; // points in which to compute polynomial
    double norm_p; // p in Lp norm
};
#endif