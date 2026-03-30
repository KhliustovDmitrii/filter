#include <cmath>

#include "polynomial_model.h"

namespace filter::examples
{
void Polynomial_Model::response(std::vector<double> &resp_arr) const
{
    int i, j;
    double xn;

    for(i=0; i<forward_size; i++)
    {
        resp_arr[i] = 0;
        xn = 1;

        for(j=num_pars-1; j>=0; j--)
        {
            resp_arr[i] = resp_arr[i] + xn*params_[j];
            xn = xn * points[i];
        }
    }
}

double Polynomial_Model::residual(std::vector<double> &mes, std::vector<double> &resp) const
{
    int i;
    double res, term;

    res = 0;
    for(i=0; i<forward_size; i++)
    {
        term = pow(std::abs(mes[i] - resp[i]), norm_p);
        res += term;
    }
    
    res = res/forward_size;
    res = pow(res, 1/norm_p);

    return res;
}

void Polynomial_Model::set_param(int ind, double val)
{
    params_[ind] = val;
}

double Polynomial_Model::get_param(int ind) const
{
    return params_[ind];
}
}; // filter::examples