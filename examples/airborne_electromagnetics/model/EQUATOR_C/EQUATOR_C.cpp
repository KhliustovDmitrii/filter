#include "EQUATOR_C.h"
#include "../../utilities/mathematics/special_functions/special_functions.h"

void EQUATOR_C::response(double *resp_arr)
{
    int i, pos;

    E->forward();

    // difference of in-phase components
    for(i=0; i<E->num_freqs-1; i++)
        resp_arr[i] = E->freq_re[i+1] - E->freq_re[i];
    
    // for highest freq duplicate the previous value
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

double EQUATOR_C::get_param(int ind)
{
    return params_[ind];
}

void EQUATOR_C::set_param(int ind, double val)
{
    double val_transformed;

    // determine which parameter it is
    std::vector<double> *dptr;
    std::vector<int> *sptr;

    int b1 = diff_rho.size();
    int b2 = b1 + diff_alt_cor.size();
    int b3 = b2 + diff_depth.size();
    int b4 = b3 + diff_cole_rho.size();
    int b5 = b4 + diff_cole_tau.size();
    int b6 = b5 + diff_cole_c.size();
    int ind_c;

    if(ind < b1) {dptr = &E->depths; sptr = &scale_rho; ind_c = ind - b1;}
    else 
    if(ind < b2) {dptr = &E->altitude_correction; sptr = &scale_alt_cor; ind_c = ind - b2;}
    else
    if(ind < b3) {dptr = &E->rhos; sptr = &scale_depth; ind_c = ind - b3;}
    else
    if(ind < b4) {dptr = &E->cole_rho; sptr = &scale_cole_rho; ind_c = ind - b4;}
    else
    if(ind < b5) {dptr = &E->cole_tau; sptr = &scale_cole_tau; ind_c = ind - b5;}
    else
    if(ind < b6) {dptr = &E->cole_c; sptr = &scale_cole_c; ind_c = ind - b6;}
    else return; // out of bounds
    
    params_[ind] = val;

    // inverse transformation given by scale
    val_transformed = inverse_transform_function_dispatch(val, sptr->at(ind_c));

    // set full model parameter
    dptr->at(ind_c) = val_transformed;
}

void EQUATOR_C::c_to_full()
{
    int i;
    double val_transformed;

    int b1 = diff_rho.size();
    int b2 = b1 + diff_alt_cor.size();
    int b3 = b2 + diff_depth.size();
    int b4 = b3 + diff_cole_rho.size();
    int b5 = b4 + diff_cole_tau.size();
    int b6 = b5 + diff_cole_c.size();
    int ind_c;

    // decompress rho
    for(i=0; i<b1; i++)
    {
        ind_c = i;
        val_transformed = inverse_transform_function_dispatch(params_[i], scale_rho[ind_c]);
        E->rhos[diff_rho[ind_c]] = val_transformed;
    }

    // decompress altitude compression
    for(; i<b2; i++)
    {
        ind_c = i-b1;
        val_transformed = inverse_transform_function_dispatch(params_[i], scale_alt_cor[ind_c]);
        E->altitude_correction[diff_alt_cor[ind_c]] = val_transformed;
    }

    // decompress depth
    for(; i<b3; i++)
    {
        ind_c = i-b2;
        val_transformed = inverse_transform_function_dispatch(params_[i], scale_depth[ind_c]);
        E->depths[diff_depth[ind_c]] = val_transformed;
    }

    // decompress cole rho
    for(; i<b4; i++)
    {
        ind_c = i-b3;
        val_transformed = inverse_transform_function_dispatch(params_[i], scale_cole_rho[ind_c]);
        E->cole_rho[diff_cole_rho[ind_c]] = val_transformed;
    }

    // decompress cole tau
    for(; i<b5; i++)
    {
        ind_c = i-b4;
        val_transformed = inverse_transform_function_dispatch(params_[i], scale_cole_tau[ind_c]);
        E->cole_tau[diff_cole_tau[ind_c]] = val_transformed;
    }

    // decompress cole c
    for(; i<b6; i++)
    {
        ind_c = i-b5;
        val_transformed = inverse_transform_function_dispatch(params_[i], scale_cole_c[ind_c]);
        E->cole_c[diff_cole_c[ind_c]] = val_transformed;
    }
}

void EQUATOR_C::full_to_c()
{
    int i;
    double val_transformed;

    int b1 = diff_rho.size();
    int b2 = b1 + diff_alt_cor.size();
    int b3 = b2 + diff_depth.size();
    int b4 = b3 + diff_cole_rho.size();
    int b5 = b4 + diff_cole_tau.size();

    for(i=0; i<diff_rho.size(); i++)
    {
        val_transformed = transform_function_dispatch(E->rhos[diff_rho[i]], scale_rho[i]);
        params_[i] = val_transformed; 
    }

    for(i=0; i<diff_alt_cor.size(); i++)
    {
        val_transformed = transform_function_dispatch(E->altitude_correction[diff_alt_cor[i]], scale_alt_cor[i]);
        params_[b1+i] = val_transformed; 
    }

    for(i=0; i<diff_depth.size(); i++)
    {
        val_transformed = transform_function_dispatch(E->depths[diff_depth[i]], scale_depth[i]);
        params_[b2+i] = val_transformed; 
    }

    for(i=0; i<diff_cole_rho.size(); i++)
    {
        val_transformed = transform_function_dispatch(E->cole_rho[diff_cole_rho[i]], scale_cole_rho[i]);
        params_[b3+i] = val_transformed; 
    }

    for(i=0; i<diff_cole_tau.size(); i++)
    {
        val_transformed = transform_function_dispatch(E->cole_tau[diff_cole_tau[i]], scale_cole_tau[i]);
        params_[b4+i] = val_transformed; 
    }

    for(i=0; i<diff_cole_c.size(); i++)
    {
        val_transformed = transform_function_dispatch(E->cole_c[diff_cole_c[i]], scale_cole_c[i]);
        params_[b5+i] = val_transformed; 
    }
}

double EQUATOR_C::residual(double *mes, double *resp)
{
    double res;
    int i;

    res = 0;
    for(i=0; i<forward_size; i++)
        res += (exp_lin(resp[i]) - exp_lin(mes[i]))*(exp_lin(resp[i]) - exp_lin(mes[i]))/(weights[i]*weights[i]);
    res = res/(forward_size);

    return sqrt(res);
}