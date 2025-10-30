#include "EQUATOR_C.h"

void aem_converter::response(double *resp_arr)
{
    int i, pos;

    E->forward();

    // difference of in-phase components
    for(i=0; i<E->num_freqs-1; i++)
        resp_arr[i] = E->freq_re[i+1] - E->freq_re[i];
    
    if(E->num_freqs > 1)
    {
        resp_arr[i] = resp_arr[i-1];
        i++;
    } 

    // quadratic components
    pos = i;
    for(i=0; i<E->num_freqs; i++)
        resp_arr[pos + i] = E->freq_im[i];

    // td channels
    pos = pos + i;
    for(i=0; i<E->num_channels; i++)
        resp_arr[pos+i] = E->td_chan[i];

    // transform response
    for(i=0; i<forward_size; i++)
        resp_arr[i] = log_lin(resp_arr[i]);
}

double get_param(int ind)
{
    return params[ind];
}

void set_param(int ind, double val)
{
    double val_transformed;
    params[ind] = val;

    switch(scale_total[i])
    {
        case 0:
            val_transformed = val;
            break;
        case 1:
            val_transformed = exp(val);
            break;
        case 2:
            val_transformed = exp_lin(val);
            break;   
    }

    if(ind < num_diff_rho)
    {
        E->rhos[ind] = val_transformed;
        return;
    }

    if(ind < num_diff_rho + diff_alt_cor)
    {
        E->altitude_correction = val_transformed;
        return;
    }

    if(ind < num_diff_rho + diff_alt_cor + num_diff_depth)
    {
        E->depths[i - num_diff_rho - diff_alt_cor] = val_transformed;
        return;
    }

    E->cole[i - num_diff_rho - diff_alt_cor - num_diff_depth] = val_transformed;
}

double update(double *upd_vec)
{
    int i;
    double val;

    for(i=0; i<num_pars; i++)
    {
        val = E->get_param(i);
        val -= upd_vec[i];
        E->set_param(i, val);
    }
}

void compress()
{
    int i;
    double val val_transformed;

    for(i=0; i<num_pars; i++)
    {
        if(ind < num_diff_rho)
            val = E->rhos[ind];
        else if(ind < num_diff_rho + diff_alt_cor)
            val = E->altitude_correction;
        else if(ind < num_diff_rho + diff_alt_cor + num_diff_depth)
            val = E->depths[i - num_diff_rho - diff_alt_cor];
        else
            val = E->cole[i - num_diff_rho - diff_alt_cor - num_diff_depth];
        
        switch(scale_total[i])
        {
            case 0:
                val_transformed = val;
                break;
            case 1:
                val_transformed = log(val);
                break;
            case 2:
                val_transformed = log_lin(val);
                break;   
        }

        params[i] = val_transformed;
    }
}


double residual(double *mes, double *resp)
{
    double res;
    int i;

    res = 0;
    for(i=0; i<forward_size; i++)
        res += (exp_lin(resp[i]) - exp_lin(mes[i]))*(exp_lin(resp[i]) - exp_lin(mes[i]))/(weights[i]*weights[i]);
    res = res/(forward_size);

    return sqrt(res);
}