#ifndef _AEM_CONVERTER_H_
#define _AEM_CONVERTER_H_

#include "../models/aem_model_raw.h"
#include "../models/model.h"
#include "../math_utils.h"

class aem_converter : public model<double, double>
{
public: 
    aem_model_raw *model;

    // differentiability - indexes of variable parameters
    int *diff_depth;
    int *diff_rho;
    int *diff_cole;
    int diff_alt_cor;

    // scales; 0 for lin, 1 for log, 2 for log-lin
    int *scale_depth;
    int *scale_rho;
    int *scale_cole;
    int scale_alt_cor;
    int *scale_total;

    // auxillary
    int num_diff_depth;
    int num_diff_rho;
    int num_diff_cole;

    double *upper_boundary;
    double *lower_boundary;

    // forward function
    void response(double *resp_arr);

    // conversion between full and filterable representation
    void raw_to_conv();
    void conv_to_raw();

    // model update method
    double update(double *upd_vec, double *resp, double res, double *R);

    aem_converter(aem_model_raw *m, int* diff_d, int* scale_d, int n_diff_d,
                                    int* diff_r, int* scale_r, int n_diff_r,
                                    int* diff_c, int* scale_c, int n_diff_c,
                                    int diff_ac, int scale_ac,
                                    double *ub, double *lb)
    {
        int i;

        model = m;

        diff_depth = diff_d;
        scale_depth = scale_d;
        num_diff_depth = n_diff_d;

        diff_rho = diff_r;
        scale_rho = scale_r;
        num_diff_rho = n_diff_r;

        diff_cole = diff_c;
        scale_cole = scale_c;
        num_diff_cole = n_diff_c;

        diff_alt_cor = diff_ac;
        scale_alt_cor = scale_ac;
        
        num_pars = num_diff_depth + num_diff_rho + num_diff_cole + diff_ac;

        forward_size = 2*m->num_freqs + m->num_channels;
        params = new double[num_pars];
        scale_total = new int[num_pars];

        for(i=0; i<num_diff_rho; i++)
            scale_total[i] = scale_rho[i];


        for(i=0; i<diff_alt_cor; i++)
            scale_total[i+num_diff_rho] = scale_alt_cor;

        for(i=0; i<num_diff_depth; i++)
            scale_total[i+num_diff_rho+diff_alt_cor] = scale_depth[i];

        for(i=0; i<num_diff_cole; i++)
            scale_total[i+num_diff_depth+num_diff_rho+diff_alt_cor] = scale_cole[i];

        upper_boundary = ub;
        lower_boundary = lb;

        raw_to_conv();
    }
};

#endif