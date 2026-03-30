#ifndef PROBABILISTIC_FILTER_H
#define PROBABILISTIC_FILTER_H

#include <vector>
#include "../model/model.h"

namespace filter
{
// computes "optimal" parameters update based on measurements
class Probabilistic_Filter
{

protected:
    // forward function computing
    Model &m;

    // noise in measurements
    std::vector<double> R;

    // noise in parameters
    std::vector<double> S;

public:

    // get updates for parameter vector and covariance matrix
    virtual void get_update(std::vector<double> &mes, 
                            std::vector<double> &upd_vec, 
                            std::vector<double> &upd_cov) const = 0;

    Probabilistic_Filter(Model &model) 
    : m(model), R(m.forward_size*m.forward_size, 0), S(m.num_pars*m.num_pars, 0) {};

    virtual ~Probabilistic_Filter() = default;

    // getters and setters for filter parameters
    std::vector<double> get_S() const;
    std::vector<double> get_R() const;

    void set_S(const std::vector<double> &S_new);
    void set_R(const std::vector<double> &R_new);
};
}; // filter
#endif