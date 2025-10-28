// g++ sensitivity.cpp ../filters/kalman_unscented.cpp ../models/aem_model_raw.cpp ../converters/aem_converter.cpp ../fft23.cpp ../text_utils.cpp ../math_utils.cpp -o forward

#include <iostream>
#include <fstream>

#include "../models/aem_model_raw.h"
#include "../converters/aem_converter.h"
#include "../filters/kalman_unscented.h"
#include "../fft23.h"
#include "../text_utils.h"
#include "../math_utils.h"


double residual(double *, double *, double *, int);

int main(int argc, char *argv[])
{
   int i, j;

   std::ifstream wf("waveform.XYZ");
   std::ofstream out(argv[1]);

   std::string buf;

   int num_freqs = 10;
   double *freqs_in = new double[num_freqs];

   double chans_in[15] = {1, 2, 3, 5, 7, 12, 19, 30, 49, 79, 128, 207, 336, 545, 884};
   int num_pol_layers = 0;
   int *pol_lays = new int[num_pol_layers+1];
   double bfr = 77.16;
   double *waveform = new double[2592];

   int num_layers = 10;

   pol_lays[0] = 0;

   for(i=0; i<num_freqs; i++) freqs_in[i] = (i+1)*100;

   aem_model_raw model = aem_model_raw(num_layers, num_pol_layers, num_freqs, 14, 100, 2592, 199, freqs_in, chans_in, pol_lays);

   for(i=1;i<2*model.num_freqs_fulltime;i+=2)
      model.freqs_fd_fulltime[(i-1)/2] = bfr*i;

    for(i=0; i<num_layers-1; i++)
        model.depths[i] = (i+1)*10;

    for(i=0; i<num_layers; i++)
        model.rhos[i] = (i+1)*100;

    model.altitude_correction = 0;

    for(i=0; i<num_pol_layers; i++)
    {
        model.cole[3*i] = 50;
        model.cole[3*i + 1] = 0.001;
        model.cole[3*i + 2] = 0.5;
    }

   model.hor_dist = 25;
   model.ver_dist = 25;
   model.alt = 50;

   model.prim_field = 1./5050000000;

   for(i=0; i<2592; i++) wf >> waveform[i];

   FFT imp_fft;
   memset(&imp_fft, 0, sizeof(imp_fft));
    
   init_fft(&imp_fft,2592);
   for(i=0;i<2592;i++)
      imp_fft.xn[i] = waveform[i];
        
   fft_pro(&imp_fft,0);
    
   for(i=0;i<2592;i++)
      model.impulse_spec[i] = imp_fft.fn[i];

   model.print_model();

   int n_diff_d = 0;
   int* diff_d = new int[n_diff_d];
   int* scale_d = new int[n_diff_d];

   int n_diff_r = num_layers;
   int *diff_r = new int[n_diff_r];
   int *scale_r = new int[n_diff_r];    
   for(i=0; i<n_diff_r; i++)       
   {
      diff_r[i] = i;
      scale_r[i] = 2;
   }

   int n_diff_c = 3;
   int *diff_c = new int[n_diff_c];
   int *scale_c = new int[n_diff_c];    
   for(i=0; i<n_diff_c; i++)       
   {
      diff_c[i] = i;
      scale_c[i] = 2;
   }

   int diff_ac = 1;
   int scale_ac = 0;

   double ub[n_diff_d + n_diff_r + n_diff_c + diff_ac];
   double lb[n_diff_d + n_diff_r + n_diff_c + diff_ac];

   for(i=0; i<n_diff_r; i++)
   {
      ub[i] = log_lin(20000);
      lb[i] = log_lin(0.01);
      scale_r[i] = 2;
      diff_r[i] = i;
   }  

   for(i=n_diff_r; i<n_diff_r+diff_ac; i++)
   {
      ub[i] = 2*model.depths[0];
      lb[i] =-2*model.depths[0];
      scale_ac = 0;
   }

   for(i=n_diff_r+diff_ac; i<n_diff_r+diff_ac+n_diff_d; i++)
   {
      ub[i] = 100;
      lb[i] = 1;
      scale_d[i - n_diff_r - diff_ac] = 0;
   } 
   
   // temporary
   int max_count = n_diff_c/3;
   for(i=0; i<max_count; i++)
   {
      ub[n_diff_d+diff_ac+n_diff_r + 3*i] = log_lin(20000);
      ub[n_diff_d+diff_ac+n_diff_r + 3*i+1] = log_lin(1);
      ub[n_diff_d+diff_ac+n_diff_r + 3*i+2] = log_lin(1);
      lb[n_diff_d+diff_ac+n_diff_r + 3*i] = log_lin(0.01);
      lb[n_diff_d+diff_ac+n_diff_r + 3*i+1] = log_lin(0.0000001);
      lb[n_diff_d+diff_ac+n_diff_r + 3*i+2] = log_lin(0.01);
      scale_c[3*i] = 2;
      scale_c[3*i+1] = 2;
      scale_c[3*i+2] = 2;
   }

   aem_converter conv = aem_converter(&model, 
                                      diff_d, scale_d, n_diff_d,
                                      diff_r, scale_r, n_diff_r,
                                      diff_c, scale_c, n_diff_c,
                                      diff_ac, scale_ac,
                                      ub, lb);

   double resp_0[conv.forward_size];
   double resp_var[conv.forward_size];
    double Hess[conv.num_pars][conv.num_pars];

    double R[conv.forward_size*conv.forward_size];
    for(i=0; i<conv.forward_size; i++)
    {
      for(j=0; j<conv.forward_size; j++)
      {
         if(i==j)
            R[conv.forward_size*i+j] = 0.1;
         else
            R[conv.forward_size*i+j] = 0;
      }
    }

    for(j=2*model.num_freqs; j<2*model.num_freqs+model.num_channels; j++)
      R[j*conv.forward_size+j] = 50000;


    // response of true model
    conv.response(resp_0);

    // residual with parameters perturbed
    double res_pp, res_pm, res_mp, res_mm;

    // with true parameters equal to 0
    double res = 0;

    double h = 0.01;

    for(i=0; i<conv.num_pars; i++)
    {
        for(j=0; j<conv.num_pars; j++)
        {
            Hess[i][j] = 0;

            conv.params[i] += h;
            conv.params[j] += h;
            conv.conv_to_raw();
            conv.response(resp_var);
            res_pp = residual(resp_var, resp_0, R, conv.forward_size);
            conv.params[i] -= h;
            conv.params[j] -= h;
            conv.conv_to_raw();

            conv.params[i] += h;
            conv.params[j] -= h;
            conv.conv_to_raw();
            conv.response(resp_var);
            res_pm = residual(resp_var, resp_0, R, conv.forward_size);
            conv.params[i] -= h;
            conv.params[j] += h;
            conv.conv_to_raw();

            conv.params[i] -= h;
            conv.params[j] += h;
            conv.conv_to_raw();
            conv.response(resp_var);
            res_mp = residual(resp_var, resp_0, R, conv.forward_size);
            conv.params[i] += h;
            conv.params[j] -= h;
            conv.conv_to_raw();

            conv.params[i] -= h;
            conv.params[j] -= h;
            conv.conv_to_raw();
            conv.response(resp_var);
            res_mm = residual(resp_var, resp_0, R, conv.forward_size);
            conv.params[i] += h;
            conv.params[j] += h;
            conv.conv_to_raw();

            Hess[i][j] = (res_pp - res_pm - res_mp - res_mm)/(4*h*h);
        }
    }

    for(i=0; i<conv.num_pars; i++)
    {
        for(j=0; j<conv.num_pars; j++)
        {
            out << Hess[i][j] << " ";
        }
        out << "\n";
    }

   return 0;
}

double residual(double *resp_1, double *resp_0, double *R, int forward_size)
{
    double res;
    int j;

    res = 0;
    for(j=0; j<forward_size; j++)
        res += (exp_lin(resp_1[j]) - exp_lin(resp_0[j]))*(exp_lin(resp_1[j]) - exp_lin(resp_0[j]))/(R[forward_size*j + j]*R[forward_size*j + j]);
    res = res/(forward_size);

    return res;
}