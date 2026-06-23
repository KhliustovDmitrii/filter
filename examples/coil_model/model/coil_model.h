#ifndef COIL_MODEL_H
#define COIL_MODEL_H

#include <vector>
#include "types/model/model.h"

namespace filter::examples
{
// model computes values of polynomial in specified ponts
// residual is calculated in Lp norm
class Coil_Model : public Model
{
public:
    Coil_Model(std::vector<double> coil_pos_1_, std::vector<double> coil_pos_2_,
    std::vector<double> mom_1_, std::vector<double> mom_2_, std::vector<double> weights_) : 
    Model(6, 6), coil_pos_1(coil_pos_1_), coil_pos_2(coil_pos_2_),
    mom_1(mom_1_), mom_2(mom_2_), params(std::vector<double>(6, 0)),
    weights(weights_) {};

    // forward function
    virtual void response(std::vector<double> &resp_arr) const override;

    // proximity in response between model and measurements
    virtual double residual(std::vector<double> &mes, std::vector<double> &resp) const override;

    // getters and setters
    virtual void set_param(int ind, double val) override;
    virtual double get_param(int ind) const override;

private:

    std::vector<double> params;

    std::vector<double> coil_pos_1; // first coil position
    std::vector<double> coil_pos_2; // second coil position
   
    std::vector<double> mom_1;      // first coil moment
    std::vector<double> mom_2;      // second coil moment

    std::vector<double> weights;    // measurement weights

    void quat_normalize(); // normalize quaternion to 1 norm
};
}; // filter::examples
#endif