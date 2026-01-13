#include "probabilistic_filter.h"

std::vector<double> Probabilistic_Filter::get_S()
{
    return S;
}

std::vector<double> Probabilistic_Filter::get_R()
{
    return R;
}

void Probabilistic_Filter::set_S(std::vector<double> S_new)
{
    S = S_new;
}

void Probabilistic_Filter::set_R(std::vector<double> R_new)
{
    R = R_new;
}