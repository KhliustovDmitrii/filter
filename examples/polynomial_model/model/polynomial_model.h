#ifndef POLYNOMIAL_MODEL_H
#define POLYNOMIAL_MODEL_H

#include <vector>
#include "types/model/model.h"

namespace filter::examples
{
// model computes values of polynomial in specified ponts
// residual is calculated in Lp norm
class Polynomial_Model : public Model
{
public:
    Polynomial_Model(std::vector<double> pts, std::vector<double> coeff, double norm) : 
    Model(coeff.size(), pts.size()), points(pts), norm_p(norm)
    {

        params = coeff;
    };

    // forward function
    virtual void response(std::vector<double> &resp_arr) const override;

    // proximity in response between model and measurements
    virtual double residual(std::vector<double> &mes, std::vector<double> &resp) const override;

    // getters and setters
    virtual void set_param(int ind, double val) override;
    virtual double get_param(int ind) const override;

private:
    std::vector<double> points; // points in which to compute polynomial
    double norm_p; // p in Lp norm

    std::vector<double> params;
};
}; // filter::examples
#endif