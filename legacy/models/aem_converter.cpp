#include "aem_converter.h"

void aem_converter::response(double *resp_arr)
{
    int i;

    model->forward();

    for(i=0; i<model->num_freqs-1; i++)
        resp_arr[2*i+1] = model->freq_re[i+1] - model->freq_re[i];
    
    if(model->num_freqs > 1) resp_arr[2*i+1] = resp_arr[2*(i-1)+1]; 

    for(i=0; i<model->num_freqs; i++)
        resp_arr[2*i] = model->freq_im[i];

    for(i=0; i<model->num_channels; i++)
        resp_arr[2*model->num_freqs+i] = model->td_chan[i];

    for(i=0; i<forward_size; i++)
        resp_arr[i] = log_lin(resp_arr[i]);
}

void aem_converter::raw_to_conv()
{
    int i;

    for(i=0; i<num_diff_rho; i++)
    {
        switch(scale_rho[i])
        {
            case 0: // linear
                params[i] = model->rhos[diff_rho[i]];
                break;
            case 1: // logarithmic
                params[i] = log(abs(model->rhos[diff_rho[i]]));
                break;
            case 2: // log-linear
                params[i] = log_lin(model->rhos[diff_rho[i]]);
                break;
        }
    }

    for(i=0; i<diff_alt_cor; i++)
    {
        switch(scale_alt_cor)
        {
            case 0: // linear
                params[i+num_diff_rho] = model->altitude_correction;
                break;
            case 1: // logarithmic
                params[i+num_diff_rho] = log(abs(model->altitude_correction));
                break;
            case 2: // log-linear
                params[i+num_diff_rho] = log_lin(model->altitude_correction);
                break;
        }
    }

    for(i=0; i<num_diff_depth; i++)
    {
        switch(scale_depth[i])
        {
            case 0: // linear
                params[i+num_diff_rho+diff_alt_cor] = model->depths[diff_depth[i]];
                break;
            case 1: // logarithmic
                params[i+num_diff_rho+diff_alt_cor] = log(abs(model->depths[diff_depth[i]]));
                break;
            case 2: // log-linear
                params[i+num_diff_rho+diff_alt_cor] = log_lin(model->depths[diff_depth[i]]);
                break;
        }
    }

    for(i=0; i<num_diff_cole; i++)
    {
        switch(scale_cole[i])
        {
            case 0: // linear
                params[i+num_diff_depth+num_diff_rho+diff_alt_cor] = model->cole[diff_cole[i]];
                break;
            case 1: // logarithmic
                params[i+num_diff_depth+num_diff_rho+diff_alt_cor] = log(abs(model->cole[diff_cole[i]]));
                break;
            case 2: // log-linear
                params[i+num_diff_depth+num_diff_rho+diff_alt_cor] = log_lin(model->cole[diff_cole[i]]);
                break;
        }
    }
}

void aem_converter::conv_to_raw()
{
    int i;

    for(i=0; i<num_diff_rho; i++)
    {
        switch(scale_rho[i])
        {
            case 0: // linear
                model->rhos[diff_rho[i]] = params[i];
                break;
            case 1: // logarithmic
                model->rhos[diff_rho[i]] = exp(params[i]);
                break;
            case 2: // log-linear
                model->rhos[diff_rho[i]] = exp_lin(params[i]);
                break;
        }
    }

    for(i=0; i<diff_alt_cor; i++)
    {
        switch(scale_alt_cor)
        {
            case 0: // linear
                model->altitude_correction = params[i+num_diff_rho];
                break;
            case 1: // logarithmic
                model->altitude_correction = exp(params[i+num_diff_rho]);
                break;
            case 2: // log-linear
                model->altitude_correction = exp_lin(params[i+num_diff_rho]);
                break;
        }
    }

    for(i=0; i<num_diff_depth; i++)
    {
        switch(scale_depth[i])
        {
            case 0: // linear
                model->depths[diff_depth[i]] = params[i+num_diff_rho+diff_alt_cor];
                break;
            case 1: // logarithmic
                model->depths[diff_depth[i]] = exp(params[i+num_diff_rho+diff_alt_cor]);
                break;
            case 2: // log-linear
                model->depths[diff_depth[i]] = exp_lin(params[i+num_diff_rho+diff_alt_cor]);
                break;
        }
    }

    for(i=0; i<num_diff_cole; i++)
    {
        switch(scale_cole[i])
        {
            case 0: // linear
                model->cole[diff_cole[i]] = params[i+num_diff_depth+num_diff_rho+diff_alt_cor];
                break;
            case 1: // logarithmic
                model->cole[diff_cole[i]] = exp(params[i+num_diff_depth+num_diff_rho+diff_alt_cor]);
                break;
            case 2: // log-linear
                model->cole[diff_cole[i]] = exp_lin(params[i+num_diff_depth+num_diff_rho+diff_alt_cor]);
                break;
        }
    }
}

double aem_converter::update(double *upd_vec, double *resp, double res, double *R)
{
    int i, count;
    double resp_tmp[forward_size];
    double old_pars[num_pars];
    double res_new;
    double factor = 1;
   
    memcpy(old_pars, params, sizeof(double)*num_pars);
    response(resp_tmp);
    res = 0;
    for(i=0; i<forward_size; i++)
        res+=(resp[i]-resp_tmp[i])*(resp[i]-resp_tmp[i])/R[i*forward_size + i];
   
    // try to reduce the residual
    for(count=0; count<3; count++)
    {
        // update model params
        for(i=0; i<num_pars; i++)
        {

            // NB! Do we need any transforms here?
            params[i] = old_pars[i] - upd_vec[i];




            
            if(std::isnan(params[i])) params[i] = 0;
            if(params[i] > upper_boundary[i]) params[i] = upper_boundary[i];
            if(params[i] < lower_boundary[i]) params[i] = lower_boundary[i];
        }
      
        conv_to_raw();

        // calculate new residual
        response(resp_tmp);
        res_new = 0;
        for(i=0; i<forward_size; i++)
            res_new+=(resp[i]-resp_tmp[i])*(resp[i]-resp_tmp[i])/R[i*forward_size + i];
      
        // try to reduce the step
        if(res_new > res)
        {
            for(i=0; i<num_pars; i++)
                upd_vec[i] = upd_vec[i]*0.5;

            factor *= 0.5;
        }
        else
        {
            count = 7;
        }
   }

   return factor;
}