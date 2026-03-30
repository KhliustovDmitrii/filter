#include "probabilistic_filter.h"

namespace filter
{
std::vector<double> Probabilistic_Filter::get_S() const
{
    return S;
}

std::vector<double> Probabilistic_Filter::get_R() const
{
    return R;
}

void Probabilistic_Filter::set_S(const std::vector<double> &S_new)
{
    S = S_new;
}

void Probabilistic_Filter::set_R(const std::vector<double> &R_new)
{
    R = R_new;
}
}; // filter