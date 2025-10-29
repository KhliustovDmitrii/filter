// g++ -g combined.cpp ../filters/kalman_filter.cpp ../filters/kalman_unscented.cpp ../models/aem_model_raw.cpp ../converters/aem_converter.cpp ../fft23.cpp ../text_utils.cpp ../math_utils.cpp -o combined

#include <iostream>
#include <fstream>

#include "../models/aem_model_raw.h"
#include "../converters/aem_converter.h"
#include "../filters/kalman_filter.h"
#include "../filters/kalman_unscented.h"
#include "../fft23.h"
#include "../text_utils.h"

#define STOP_VAL 1.
#define rho_0 10

void set_model(aem_model_raw *);
double residual(double *, double *, double *, int);

int main()
{
   std::ifstream data_file("Italy_short.XYZ");
   std::ofstream output("10_ext_3_uns.XYZ");
   std::string buf;
   int i;
   // -------------------------------------- 1. Set up models
   int num_layers = 25;
   int num_pol_layers = 0;
   int num_freqs = 21;
   int num_channels = 14;
   double bfr = 77.16;
   double err_ini = 0.3;
   double cor_ini = 0.1;

   double freqs_in[num_freqs] = {77.1, 231.5, 385.8, 540.1, 694.4, 1003.1, 1157.41, 1311.7, 1466, 1620.4, 1774.7,
                       1929, 2083.3, 2237.65, 2854.94, 3009.26, 3163.58, 3317.9, 3472.22, 6172.84, 13503.09};
   double chans_in[15] = {1, 2, 3, 5, 7, 12, 19, 30, 49, 79, 128, 207, 336, 545, 884};
   int pol_layers[num_pol_layers+1];
   double noise_fd[num_freqs] = {0.2, 0.12, 0.09, 0.13, 0.39, 0.5, 0.58, 0.52, 0.36, 0.46, 0.15, 0.31, 1.01, 0.56, 0.63,
                       0.76, 0.44, 0.23, 0.44, 0.52, 1.11};
   double noise_td[num_channels] = {71000, 71000, 70000, 66000, 53000, 33000, 33000, 23000, 19000, 30000, 18000, 16000, 10000, 7000};
   
   aem_model_raw model_f = aem_model_raw(num_layers, num_pol_layers, num_freqs, num_channels,
                                       100, 2592, 199, freqs_in, chans_in, pol_layers);
   set_model(&model_f);

   double freqs_red[2] = {77.1, 231.5};
   int empty_pol_layers[1];
   aem_model_raw model_r = aem_model_raw(1, 0, 2, 0, 
                                        100, 2592, 199, freqs_red, chans_in, empty_pol_layers);
   model_r.altitude_correction = 0;
   model_r.prim_field = primField(32, 25)/10000;

   
   // -------------------------------------- 2. Set up converter for full model
   int n_diff_d = 0;
   int diff_d[n_diff_d];
   int scale_d[n_diff_d];

   int n_diff_r = 25;
   int diff_r[n_diff_r];
   int scale_r[n_diff_r];                               

   int n_diff_c = 0;
   int diff_c[3*n_diff_c];
   int scale_c[3*n_diff_c];

   int diff_ac = 1;
   int scale_ac = 0;

   double ub[n_diff_d + n_diff_r + n_diff_c + diff_ac];
   double lb[n_diff_d + n_diff_r + n_diff_c + diff_ac];

   for(i=0; i<n_diff_r; i++)
   {
      ub[i] = log_lin(20000);
      lb[i] = log_lin(0.001);
      scale_r[i] = 2;
      diff_r[i] = i;
   } 

   for(i=n_diff_r; i<n_diff_r+diff_ac; i++)
   {
      ub[i] = 2*model_f.depths[0];
      lb[i] =-2*model_f.depths[0];
      scale_ac = 0;
   }

   aem_converter conv_f = aem_converter(&model_f, 
                                      diff_d, scale_d, n_diff_d,
                                      diff_r, scale_r, n_diff_r,
                                      diff_c, scale_c, n_diff_c,
                                      diff_ac, scale_ac,
                                      ub, lb);

   // -------------------------------------- 3. Set up converter for reduced model
   int n_diff_d_s = 0;
   int diff_d_s[n_diff_d_s];
   int scale_d_s[n_diff_d_s];

   int n_diff_r_s = 1;
   int diff_r_s[n_diff_r_s];
   int scale_r_s[n_diff_r_s];                               

   int n_diff_c_s = 0;
   int diff_c_s[3*n_diff_c_s];
   int scale_c_s[3*n_diff_c_s];

   int diff_ac_s = 0;
   int scale_ac_s = 0;

   double ub_s[n_diff_d_s + n_diff_r_s + n_diff_c_s + diff_ac_s];
   double lb_s[n_diff_d_s + n_diff_r_s + n_diff_c_s + diff_ac_s];

   for(i=0; i<n_diff_r_s; i++)
   {
      ub_s[i] = log_lin(20000);
      lb_s[i] = log_lin(0.001);
      diff_r_s[i] = i;
      scale_r_s[i] = 2;
   }  

   aem_converter conv_r = aem_converter(&model_r, 
                                      diff_d_s, scale_d_s, n_diff_d_s,
                                      diff_r_s, scale_r_s, n_diff_r_s,
                                      diff_c_s, scale_c_s, n_diff_c_s,
                                      diff_ac_s, scale_ac_s,
                                      ub_s, lb_s);

   // -------------------------------------- 3. Set up filters
   kalman_filter k_r = kalman_filter(&conv_r);
   kalman_filter k_f_ext = kalman_filter(&conv_f);
   kalman_unscented k_f_uns = kalman_unscented(&conv_f);

   int forward_size = k_f_ext.m->forward_size;
   int num_pars = k_f_ext.m->num_pars;

   k_r.R = new double[k_r.m->forward_size*k_r.m->forward_size];
   k_r.S = new double[k_r.m->num_pars*k_r.m->num_pars];

   k_f_ext.R = new double[forward_size*forward_size];
   k_f_ext.S = new double[num_pars*num_pars];
   k_f_uns.R = new double[forward_size*forward_size];
   k_f_uns.S = new double[num_pars*num_pars];

   double R_0[forward_size*forward_size];
   memset(R_0, 0, forward_size*forward_size*sizeof(double));
   for(i=0; i<model_f.num_freqs; i++)
      R_0[forward_size*2*i+2*i] = noise_fd[i];

   for(i=0; i<model_f.num_freqs-1; i++)
      R_0[forward_size*(2*i+1)+2*i+1] = sqrt(noise_fd[i]*noise_fd[i] + noise_fd[i+1]*noise_fd[i+1]);

   R_0[forward_size*(2*i+1)+2*i+1] = 1000;

   for(i=2*model_f.num_freqs; i<2*model_f.num_freqs+model_f.num_channels; i++)
      R_0[forward_size*i+i] = noise_td[i-2*model_f.num_freqs];

   double S_0[num_pars*num_pars];
   for(i=num_layers-1;i>=0;i--) 
   {
      if(i==num_layers-1)//Bottom layer resisitivity
         S_0[i+num_pars*i] = err_ini;
      else //Other resisitivities
      {
         S_0[i+1+num_pars*i] = cor_ini*err_ini/S_0[i+1+num_pars*(i+1)];
         S_0[i+num_pars*i] = sqrt(err_ini*err_ini-S_0[i+1+num_pars*i]*S_0[i+1+num_pars*i]);
      }
   }

   // ac 
   S_0[num_layers + num_layers*num_pars] = err_ini;

   double R_0_r[k_r.m->forward_size*k_r.m->forward_size];
   memset(R_0_r, 0, k_r.m->forward_size*k_r.m->forward_size*sizeof(double));
   R_0_r[0] = R_0[0];
   R_0_r[k_r.m->forward_size + 1] = R_0[forward_size + 1];
   R_0_r[2*k_r.m->forward_size + 2] = 1000;
   R_0_r[3*k_r.m->forward_size + 3] = 1000;

   double S_0_r[k_r.m->num_pars*k_r.m->num_pars];
   S_0_r[0] = 0.1;

   double resp[forward_size];
   double mes_parsed[forward_size];
   double mes[forward_size];
   double mes_reduced[k_r.m->forward_size];

   int agr_count = 0;
   double factor;
   int hd_pos = 0;
   int vd_pos = 1;
   int alt_pos = 2;
   int mes_pos = 3;
   int size_of_time = 17;
   double res_prev = 1000;
   int started = 1;
   double weight;
   std::string time;
   double summed_vals[256];
   double parsed_vals[256];
   double tmp_vals[256];
   memset(summed_vals, 0, 256*sizeof(double));
   memset(parsed_vals, 0, 256*sizeof(double));
   memset(tmp_vals, 0, 256*sizeof(double));

   double best_pars[num_pars];
   double upd_vec[num_pars];
   double upd_cov[num_pars*num_pars];
   double uncert_buf[num_pars];

   double res, res_best;

   int iter;

   while(getline(data_file, buf))
   {
      parse(&(buf[size_of_time]), tmp_vals);
      time = buf.substr(0, size_of_time);

      if(agr_count < 4) 
      {
         for(i=0; i<256; i++)
            summed_vals[i] += tmp_vals[i];

         agr_count++;
         continue;
      } 
      else
      {
         for(i=0; i<256; i++)
            parsed_vals[i] = summed_vals[i]/4;

         memset(summed_vals, 0, 256*sizeof(double));
         for(i=0; i<256; i++)
            summed_vals[i] += tmp_vals[i];

         agr_count = 1;
      }

      std::cout << time << "\n";

      memcpy(mes_parsed, &(parsed_vals[mes_pos]), forward_size*sizeof(double));

      // make re positive
      for(i=0; i<num_freqs; i++)
         mes_parsed[i] = - mes_parsed[i];

      // re to diff
      for(i=0; i<model_f.num_freqs-1; i++)
         mes[2*i+1] = mes_parsed[i+1] - mes_parsed[i];

      mes[2*i+1] = mes[2*i-1];

      // im
      for(i=0; i<model_f.num_freqs; i++)
         mes[2*i] = mes_parsed[i+model_f.num_freqs];

      // td
      for(i=2*model_f.num_freqs; i<forward_size; i++)
         mes[i] = mes_parsed[i];

      for(i=0; i<forward_size; i++) 
         mes[i] = log_lin(mes[i]);

      for(i=0; i<k_r.m->forward_size; i++)
         mes_reduced[i] = mes[i];

      mes_reduced[k_r.m->forward_size - 1] = mes_reduced[k_r.m->forward_size - 3];

      std::cout << "FITTING REDUCED MODEL...\n";

      model_r.rhos[0] = rho_0;
      model_r.hor_dist = parsed_vals[hd_pos];
      model_r.ver_dist = parsed_vals[vd_pos];
      model_r.alt = parsed_vals[alt_pos];
      conv_r.raw_to_conv();

      k_r.m->response(resp);
      res_best = residual(resp, mes_reduced, R_0_r, k_r.m->forward_size);
      memcpy(k_r.R, R_0_r, k_r.m->forward_size*k_r.m->forward_size*sizeof(double));
      for(iter=0; iter<10; iter++)
      {
         memcpy(k_r.S, S_0_r, k_r.m->num_pars * k_r.m->num_pars * sizeof(double));
         memcpy(best_pars, k_r.m->params, k_r.m->num_pars*sizeof(double));
         k_r.get_update_vector(upd_vec, mes_reduced);      
         k_r.m->update(upd_vec, mes_reduced, res, R_0_r);
         k_r.m->response(resp);
         
         res = residual(resp, mes_reduced, R_0_r, k_r.m->forward_size);
         std::cout << "-------------- " << res << "\n";

         if(res < 1) break;
         if(res > res_best) 
         {
            memcpy(k_r.m->params, best_pars, k_r.m->num_pars*sizeof(double));
            conv_r.conv_to_raw();
            res = res_best;
            break;
         }

         res_best = res;
      }

      std::cout << "REDUCED MODEL FITTED\n";
      model_r.print_model();

      std::cout << "FITTING FULL MODEL...\n";

      model_f.hor_dist = parsed_vals[hd_pos];
      model_f.ver_dist = parsed_vals[vd_pos];
      model_f.alt = parsed_vals[alt_pos];

      if(started == 1)
      {
         for(i=0; i<num_layers; i++)
            model_f.rhos[i] = model_r.rhos[0];

         started = 0;
      }
      else
      {
         weight = std::min(1./res_prev, 1.);
         for(i=0; i<num_layers; i++)
            model_f.rhos[i] = weight*model_f.rhos[i] + (1-weight)*model_r.rhos[0];
      }

      conv_f.raw_to_conv();
      k_f_ext.m->response(resp);
      res_best = residual(resp, mes, R_0, forward_size);

      // FIT EXTENDED FILTER
      // adjust for log-lin
      for(i=0; i<forward_size; i++) 
      {
         if (abs(mes[i]) > 5*R_0[i*forward_size + i])
            k_f_ext.R[i*forward_size + i] = R_0[i*forward_size + i]/mes[i];
         else
            k_f_ext.R[i*forward_size + i] = R_0[i*forward_size + i];
      }

      for(iter=0; iter<10; iter++)
      {
         memcpy(k_f_ext.S, S_0, num_pars * num_pars * sizeof(double));
         memcpy(best_pars, k_f_ext.m->params, k_f_ext.m->num_pars*sizeof(double));
         k_f_ext.get_update_vector(upd_vec, mes);      
         k_f_ext.m->update(upd_vec, mes, res, k_f_ext.R);
         k_f_ext.m->response(resp);
         
         res = residual(resp, mes, R_0, forward_size);
         std::cout << "-------------- " << res << "\n";

         if(res < 1) break;
         if(res > res_best) 
         {
            memcpy(k_f_ext.m->params, best_pars, num_pars*sizeof(double));
            conv_f.conv_to_raw();
            res = res_best;
            break;
         }

         res_best = res;
      }

      std::cout << "FULL MODEL FITTED\n";
      model_f.print_model();

      std::cout << "TUNING MODEL...\n";

      // FIT UNSCENTED FILTER
      // adjust for log-lin
      for(i=0; i<forward_size; i++) 
      {
         if (abs(mes[i]) > 5*R_0[i*forward_size + i])
            k_f_uns.R[i*forward_size + i] = R_0[i*forward_size + i]/mes[i];
         else
            k_f_uns.R[i*forward_size + i] = R_0[i*forward_size + i];
      }

      memcpy(k_f_uns.S, S_0, num_pars * num_pars * sizeof(double));
      for(iter=0; iter<3; iter++)
      {
         memcpy(best_pars, k_f_uns.m->params, k_f_uns.m->num_pars*sizeof(double));
         k_f_uns.compute_sigma();  
         k_f_uns.get_matrices();  
         k_f_uns.get_update_vector(upd_vec, upd_cov, mes);      
         factor = k_f_uns.m->update(upd_vec, mes, res, k_f_uns.R);
         k_f_uns.update_cov(upd_cov, factor);
         k_f_uns.m->response(resp);
         res = residual(resp, mes, R_0, forward_size);
         std::cout << "---- " << factor << "\n";
         std::cout << "-------------- " << res << "\n";

         if(res < 1) break;
         
         if(res > res_best) 
         {
            memcpy(k_f_uns.m->params, best_pars, num_pars*sizeof(double));
            conv_f.conv_to_raw();
            res = res_best;
            break;
         }
         
         res_best = res;
      }

      std::cout << "FULL MODEL TUNED\n";
      model_f.print_model();

      output << time << " ";
      output << res << " ";
      model_f.print_to_file(output);

      output << model_r.rhos[0] << " ";

      var_mes(S_0, k_f_uns.S, uncert_buf, num_pars);

      for(i=0; i<num_pars; i++) 
         output << uncert_buf[i] << " ";

      output << std::endl;
      res_prev = res;
   }

   return 0;
}

void set_model(aem_model_raw *model)
{
   std::ifstream wf("waveform.XYZ");
   double bfr = 77.16;
   double *waveform = new double[2592];
   int i;
   
   for(i=0; i<2592; i++) wf >> waveform[i];

   FFT imp_fft;
   memset(&imp_fft, 0, sizeof(imp_fft));
    
   init_fft(&imp_fft,2592);
   for(i=0;i<2592;i++)
      imp_fft.xn[i] = waveform[i];
        
   fft_pro(&imp_fft,0);

   for(i=1;i<2*model->num_freqs_fulltime;i+=2)
      model->freqs_fd_fulltime[(i-1)/2] = bfr*i;

   model->depths[0] = 4;
   for(i=1; i<25; i++)
      model->depths[i] = model->depths[i-1]*1.1085;

   model->altitude_correction = 0;

   model->prim_field = primField(32, 25)/10000;
    
   for(i=0;i<2592;i++)
      model->impulse_spec[i] = imp_fft.fn[i];
}

double residual(double *resp_1, double *resp_0, double *R, int forward_size)
{
    double res;
    int j;

    res = 0;
    for(j=0; j<forward_size; j++)
        res += (exp_lin(resp_1[j]) - exp_lin(resp_0[j]))*(exp_lin(resp_1[j]) - exp_lin(resp_0[j]))/(R[forward_size*j + j]*R[forward_size*j + j]);
    res = res/(forward_size);

    return sqrt(res);
}


