#ifndef _EQUATOR_C_H_
#define _EQUATOR_C_H_

// compressed representation of EQUATOR
// complies to model interface
#include "../EQUATOR/EQUATOR.h"
#include "types/model/model.h"

class EQUATOR_data_loader;
class EQUATOR_aggregator;

class EQUATOR_C : public Model
{

public:

    friend class EQUATOR_data_loader;
    friend class EQUATOR_aggregator;

    EQUATOR_C(EQUATOR *EQ,
              std::vector<int> dd, std::vector<int> dr, std::vector<int> dac,
              std::vector<int> dcr, std::vector<int> dct, std::vector<int> dcc,
              std::vector<int> sd, std::vector<int> sr, std::vector<int> sac,
              std::vector<int> scr, std::vector<int> sct, std::vector<int> scc,
              std::vector<double> w) :
              Model(),
              E(EQ), diff_depth(dd), diff_rho(dr), diff_alt_cor(dac),
              diff_cole_rho(dcr), diff_cole_tau(dct), diff_cole_c(dcc),
              scale_depth(sd), scale_rho(sr), scale_alt_cor(sac),
              scale_cole_rho(scr), scale_cole_tau(sct), scale_cole_c(scc),
              weights(w)
              {

                // components of the interface
                num_pars = diff_depth.size() + diff_rho.size() + diff_alt_cor.size() +
                diff_cole_rho.size() + diff_cole_tau.size() + diff_cole_c.size();

                forward_size = 2 * E->num_freqs + E->num_channels;

                params_.resize(num_pars, 0.0);

                // set up params equal to full EQUATOR's
                full_to_c();
              };

    // forward function
    virtual void response(std::vector<double> &resp_arr) override;

    // proximity in response between model and measurements
    virtual double residual(std::vector<double> &mes, std::vector<double> &resp) override;

    // getters and setters
    virtual void set_param(int ind, double val) override;
    virtual double get_param(int ind) override;

private:

    EQUATOR *E;

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

#endif