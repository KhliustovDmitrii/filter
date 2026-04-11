#include "EQUATOR_C.h"
#include "utilities/mathematics/special_functions/special_functions.h"

namespace filter::examples
{
void EQUATOR_C::response(std::vector<double> &resp_arr) const
{
    int i, pos;

    E.forward();

    // difference of in-phase components
    for(i=0; i<E.num_freqs-1; i++)
        resp_arr[i] = E.freq_re[i+1] - E.freq_re[i];
    
    // for highest freq duplicate the previous value
    if(E.num_freqs > 1)
    {
        resp_arr[i] = resp_arr[i-1];
        i++;
    } 

    // quadratic components
    pos = i;
    for(i=0; i<E.num_freqs; i++)
        resp_arr[pos + i] = E.freq_im[i];

    // td channels
    pos = pos + i;
    for(i=0; i<E.num_channels; i++)
        resp_arr[pos+i] = E.td_chan[i];

    // transform response
    for(i=0; i<forward_size; i++)
        resp_arr[i] = filter::math::log_lin(resp_arr[i]);
}

double EQUATOR_C::get_param(int ind) const
{
    return params[ind]->get();
}

void EQUATOR_C::set_param(int ind, double val)
{
    params[ind]->set(val);
}

double EQUATOR_C::residual(std::vector<double> &mes, std::vector<double> &resp) const
{
    double res;
    int i;

    res = 0;
    for(i=0; i<forward_size; i++)
        res += (filter::math::exp_lin(resp[i]) - filter::math::exp_lin(mes[i]))*(filter::math::exp_lin(resp[i]) - filter::math::exp_lin(mes[i]))/(weights[i]*weights[i]);
    res = res/(forward_size);

    return sqrt(res);
}

}; // filter::examples