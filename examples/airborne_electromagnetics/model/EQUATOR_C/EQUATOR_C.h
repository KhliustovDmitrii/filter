#ifndef EQUATOR_C_H
#define EQUATOR_C_H

// compressed representation of EQUATOR
// complies to model interface
#include "../EQUATOR/EQUATOR.h"
#include "types/model/model.h"

namespace filter::examples
{
class EQUATOR_data_loader;
class EQUATOR_aggregator;

class EQUATOR_C : public Model
{

public:

    friend class EQUATOR_data_loader;
    friend class EQUATOR_aggregator;

    EQUATOR_C(EQUATOR &EQ,
              std::vector<int> diff_depth_, std::vector<int> diff_rho_, 
              std::vector<int> diff_alt_cor_,
              std::vector<int> diff_cole_rho_, std::vector<int> diff_cole_tau_, 
              std::vector<int> diff_cole_c_,
              std::vector<int> scale_depth_, std::vector<int> scale_rho_, 
              std::vector<int> scale_alt_cor_,
              std::vector<int> scale_cole_rho_, std::vector<int> scale_cole_tau_, 
              std::vector<int> scale_cole_c_,
              std::vector<double> weights_) :
              Model(
                diff_depth_.size() + diff_rho_.size() + diff_alt_cor_.size() +
                diff_cole_rho_.size() + diff_cole_tau_.size() + diff_cole_c_.size(),
                2 * E.num_freqs + E.num_channels),
              E(EQ), diff_depth(diff_depth_), diff_rho(diff_rho_), diff_alt_cor(diff_alt_cor_),
              diff_cole_rho(diff_cole_rho_), diff_cole_tau(diff_cole_tau_), diff_cole_c(diff_cole_c_),
              scale_depth(scale_depth_), scale_rho(scale_rho_), scale_alt_cor(scale_alt_cor_),
              scale_cole_rho(scale_cole_rho_), scale_cole_tau(scale_cole_tau_), scale_cole_c(scale_cole_c_),
              weights(weights_)
              {
                // set up params equal to full EQUATOR's
                full_to_c();
              };

    // forward function
    virtual void response(std::vector<double> &resp_arr) const override;

    // proximity in response between model and measurements
    virtual double residual(std::vector<double> &mes, std::vector<double> &resp) const override;

    // getters and setters
    virtual void set_param(int ind, double val) override;
    virtual double get_param(int ind) const override;

private:

    EQUATOR &E;

    // differentiability of EQUATOR model parameters
    std::vector<int> diff_depth, diff_rho;
    std::vector<int> diff_cole_rho, diff_cole_tau, diff_cole_c;

    // ugly, but we use vector here for uniformity
    // empty corresponds to non-differentiable parameter
    std::vector<int> diff_alt_cor;

    // scales; 0 for lin, 1 for log, 2 for log-lin
    std::vector<int> scale_depth, scale_rho;
    std::vector<int> scale_cole_rho, scale_cole_tau, scale_cole_c;
    std::vector<int> scale_alt_cor;

    // weights to compute residual with
    std::vector<double> weights;

    // transforms underlying EQUATOR to _C  and vice versa
    void c_to_full();
    void full_to_c();
};
}; // filter::examples
#endif