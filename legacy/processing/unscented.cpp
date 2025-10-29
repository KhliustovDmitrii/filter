// g++ unscented.cpp ../filters/kalman_unscented.cpp ../models/aem_model_raw.cpp ../converters/aem_converter.cpp ../fft23.cpp ../text_utils.cpp ../math_utils.cpp -o unscented

#include <iostream>
#include <fstream>

#include "../models/aem_model_raw.h"
#include "../converters/aem_converter.h"
#include "../filters/kalman_unscented.h"
#include "../fft23.h"
#include "../text_utils.h"

#define STOP_VAL 1.
#define rho_0 10

int main()
{
   int i, j, it, l;
   double res, res_old;

   std::ifstream data_file("Italy_short.XYZ");
   std::ofstream output("result_unscent.XYZ");
   std::ifstream wf("waveform.XYZ");
   std::string buf;

   // -------------------------------------- 1. Read spectrum
   double *waveform = new double[2592];
   
   for(i=0; i<2592; i++) wf >> waveform[i];

   FFT imp_fft;
   memset(&imp_fft, 0, sizeof(imp_fft));
    
   init_fft(&imp_fft,2592);
   for(i=0;i<2592;i++)
      imp_fft.xn[i] = waveform[i];
        
   fft_pro(&imp_fft,0);

   // -------------------------------------- 2. Set up full model
   int num_layers = 25;
   int num_pol_layers = 0;
   int num_freqs = 21;
   int num_channels = 14;

   double freqs_in[num_freqs] = {77.1, 231.5, 385.8, 540.1, 694.4, 1003.1, 1157.41, 1311.7, 1466, 1620.4, 1774.7,
                       1929, 2083.3, 2237.65, 2854.94, 3009.26, 3163.58, 3317.9, 3472.22, 6172.84, 13503.09};
   double freqs_red[2] = {77.1, 231.5};
   double chans_in[15] = {1, 2, 3, 5, 7, 12, 19, 30, 49, 79, 128, 207, 336, 545, 884};
   int pol_layers[num_pol_layers+1];
   int empty_pol_layers[1];
   double bfr = 77.16;

   aem_model_raw model = aem_model_raw(num_layers, num_pol_layers, num_freqs, num_channels,
                                       100, 2592, 199, freqs_in, chans_in, pol_layers);
   aem_model_raw simple_model = aem_model_raw(1, 0, 2, 0, 
                                              100, 2592, 199, freqs_red, chans_in, empty_pol_layers);

   for(i=1;i<2*model.num_freqs_fulltime;i+=2)
      model.freqs_fd_fulltime[(i-1)/2] = bfr*i;

   model.depths[0] = 4;
   for(i=1; i<25; i++)
      model.depths[i] = model.depths[i-1]*1.1085;

   model.altitude_correction = 0;

   //model.prim_field = 1./5050000000;
   model.prim_field = primField(32, 25)/10000;
    
   for(i=0;i<2592;i++)
      model.impulse_spec[i] = imp_fft.fn[i];

   // -------------------------------------- 2. Set up reduced model
   simple_model.rhos[0] = rho_0;
   
   simple_model.altitude_correction = 0;

   //simple_model.prim_field = 1./5050000000;
   simple_model.prim_field = primField(32, 25)/10000;

   // -------------------------------------- 3. Set up converter for full model
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
      ub[i] = 20000;
      lb[i] = 0.001;
      scale_r[i] = 2;
      diff_r[i] = i;
   } 

   for(i=n_diff_r; i<n_diff_d+diff_ac; i++)
   {
      ub[i] = 2*model.depths[0];
      lb[i] =-2*model.depths[0];
      scale_ac = 0;
   }

   aem_converter conv = aem_converter(&model, 
                                      diff_d, scale_d, n_diff_d,
                                      diff_r, scale_r, n_diff_r,
                                      diff_c, scale_c, n_diff_c,
                                      diff_ac, scale_ac,
                                      ub, lb);

   // -------------------------------------- 4. Set up converter for reduced model
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
      ub_s[i] = 20000;
      lb_s[i] = 0.001;
      diff_r_s[i] = i;
      scale_r_s[i] = 2;
   }  

   aem_converter simple_conv = aem_converter(&simple_model, 
                                      diff_d_s, scale_d_s, n_diff_d_s,
                                      diff_r_s, scale_r_s, n_diff_r_s,
                                      diff_c_s, scale_c_s, n_diff_c_s,
                                      diff_ac_s, scale_ac_s,
                                      ub_s, lb_s);

   // -------------------------------------- 5. Set up filter for full model
   kalman_unscented k = kalman_unscented(&conv);

   k.R = new double[k.m->forward_size*k.m->forward_size];
   k.S = new double[k.m->num_pars*k.m->num_pars];

   for(i=0; i<k.m->num_pars; i++)
   {
      for(j=0; j<k.m->num_pars; j++)
      {
         if(i==j)
            k.S[k.m->num_pars*i+j] = 1;
         else
            k.S[k.m->num_pars*i+j] = 0;
      }
   }

   double err_ini = 0.3;
   double cor_ini = 0.1;
   for(i=model.num_layers-1;i>=0;i--) 
   {
      if(i==model.num_layers-1)//Bottom layer resisitivity
         k.S[i+k.m->num_pars*i] = err_ini;
      else //Other resisitivities
      {
         k.S[i+1+k.m->num_pars*i] = cor_ini*err_ini/k.S[i+1+k.m->num_pars*(i+1)];
         k.S[i+k.m->num_pars*i] = sqrt(err_ini*err_ini-k.S[i+1+k.m->num_pars*i]*k.S[i+1+k.m->num_pars*i]);
      }
   }

   memset(k.R, 0, k.m->forward_size*k.m->forward_size*sizeof(double));

   double noise_fd[num_freqs] = {0.2, 0.12, 0.09, 0.13, 0.39, 0.5, 0.58, 0.52, 0.36, 0.46, 0.15, 0.31, 1.01, 0.56, 0.63,
                       0.76, 0.44, 0.23, 0.44, 0.52, 1.11};
   double noise_td[num_channels] = {71000, 71000, 70000, 66000, 53000, 33000, 33000, 23000, 19000, 30000, 18000, 16000, 10000, 7000};

   for(i=0; i<model.num_freqs; i++)
      k.R[k.m->forward_size*2*i+2*i] = noise_fd[i];

   for(i=0; i<model.num_freqs-1; i++)
      k.R[k.m->forward_size*(2*i+1)+2*i+1] = sqrt(noise_fd[i]*noise_fd[i] + noise_fd[i+1]*noise_fd[i+1]);

   k.R[k.m->forward_size*(2*i+1)+2*i+1] = 1000;

   for(i=2*model.num_freqs; i<2*model.num_freqs+model.num_channels; i++)
      k.R[i*k.m->forward_size+i] = noise_td[i-2*model.num_freqs];

   // -------------------------------------- 6. Set up filter for reduced model
   kalman_unscented k_s = kalman_unscented(&simple_conv);

   k_s.R = new double[k_s.m->forward_size*k_s.m->forward_size];
   k_s.S = new double[k_s.m->num_pars*k_s.m->num_pars];

   for(i=0; i<k_s.m->num_pars; i++)
   {
      for(j=0; j<k_s.m->num_pars; j++)
      {
         if(i==j)
            k_s.S[k_s.m->num_pars*i+j] = 1;
         else
            k_s.S[k_s.m->num_pars*i+j] = 0;
      }
   }

   memset(k_s.R, 0, k_s.m->forward_size*k_s.m->forward_size*sizeof(double));

   for(i=0; i<simple_model.num_freqs; i++)
      k_s.R[k_s.m->forward_size*2*i+2*i] = noise_fd[i];

   for(i=0; i<simple_model.num_freqs-1; i++)
      k_s.R[k_s.m->forward_size*(2*i+1)+2*i+1] = sqrt(noise_fd[i]*noise_fd[i] + noise_fd[i+1]*noise_fd[i+1]);

   k_s.R[k_s.m->forward_size*(2*i+1)+2*i+1] = 1000;

   k_s.R[k_s.m->forward_size*(2*i)+2*i] = 1000;

   for(i=2*simple_model.num_freqs; i<2*simple_model.num_freqs+simple_model.num_channels; i++)
      k_s.R[i*k_s.m->forward_size+i] = noise_td[i-2*simple_model.num_freqs];

   double resp[k.m->forward_size];
   double resp_reduced[k_s.m->forward_size];
   double resp_tmp[k.m->forward_size];
   double resp_parsed[k.m->forward_size];
   double best_params[k.m->num_pars];
   double upd_vec[k.m->num_pars];
   double upd_cov[k.m->num_pars*k.m->num_pars];
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
   double uncert_buf[k.m->num_pars];
   double S0[k.m->num_pars * k.m->num_pars];
   double S0_s[k_s.m->num_pars * k_s.m->num_pars];
   double old_pars[k.m->num_pars];
   double R_adj[k.m->forward_size*k.m->forward_size];
   double R_0[k.m->forward_size*k.m->forward_size];

   // base noise of filter
   memcpy(R_0, k.R, k.m->forward_size*k.m->forward_size*sizeof(double));

   // noise adjusted for log_lin
   memcpy(R_adj, k.R, k.m->forward_size*k.m->forward_size*sizeof(double));

   memcpy(S0, k.S, k.m->num_pars * k.m->num_pars * sizeof(double));
   memcpy(S0_s, k_s.S, k_s.m->num_pars * k_s.m->num_pars * sizeof(double));
   // -------------------------------------- 4. Start the data processing cycle
   int agr_count = 0;
   double summed_vals[256];
   double parsed_vals[256];
   double tmp_vals[256];
   memset(summed_vals, 0, 256*sizeof(double));
   memset(parsed_vals, 0, 256*sizeof(double));
   memset(tmp_vals, 0, 256*sizeof(double));
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

      memcpy(k.S, S0, k.m->num_pars * k.m->num_pars * sizeof(double));
      memcpy(k_s.S, S0_s, k_s.m->num_pars * k_s.m->num_pars * sizeof(double));

      parse(&(buf[size_of_time]), parsed_vals);
      time = buf.substr(0, size_of_time);

      simple_model.hor_dist = parsed_vals[hd_pos];
      simple_model.ver_dist = parsed_vals[vd_pos];
      simple_model.alt = parsed_vals[alt_pos];

      memcpy(resp_parsed, &(parsed_vals[mes_pos]), k.m->forward_size*sizeof(double));

      memcpy(resp, resp_parsed, k.m->forward_size*sizeof(double));
      // make re positive
      for(i=0; i<num_freqs; i++)
         resp_parsed[i] = - resp_parsed[i];

      // re to diff
      for(i=0; i<model.num_freqs-1; i++)
         resp[2*i+1] = resp_parsed[i+1] - resp_parsed[i];

      resp[2*i+1] = resp[2*(i-1)+1];

      // im
      for(i=0; i<model.num_freqs; i++)
         resp[2*i] = resp_parsed[i+model.num_freqs];

      for(i=0; i<conv.forward_size; i++) 
         resp[i] = log_lin(resp[i]);

      for(i=0; i<simple_conv.forward_size; i++) 
         resp_reduced[i] = resp[i];

      if(simple_model.num_freqs > 1) resp_reduced[2*simple_model.num_freqs-1] = resp_reduced[2*simple_model.num_freqs-3];
      // -------------------------------------- 4.1. Fit the simple model

      std::cout << "FITTING REDUCED MODEL...\n";
      simple_model.rhos[0] = rho_0;

      k_s.m->response(resp_tmp);
      res_old = 0;
      for(j=0; j<k_s.m->forward_size; j++)
         res_old += (exp_lin(resp_tmp[j]) - exp_lin(resp_reduced[j]))*(exp_lin(resp_tmp[j]) - exp_lin(resp_reduced[j]))/(k_s.R[k_s.m->forward_size*j + j]*k_s.R[k_s.m->forward_size*j + j]);

      res_old = res_old/(2*simple_model.num_freqs + simple_model.num_channels);
      std::cout << "-------------+ " << res_old << "\n";

      
      memcpy(k_s.S, S0_s, k_s.m->num_pars * k_s.m->num_pars * sizeof(double));
      for(it = 0; it < 10; it++)
      {
         memcpy(old_pars, k_s.m->params, k_s.m->num_pars*sizeof(double));
         k_s.compute_sigma();  
         k_s.get_matrices();  
         k_s.get_update_vector(upd_vec, upd_cov, resp_reduced);      
         factor = k_s.m->update(upd_vec, resp_reduced, res, k_s.R);
         k_s.update_cov(upd_cov, factor);
         
         k_s.m->response(resp_tmp);
         
         res = 0;
         for(j=0; j<k_s.m->forward_size; j++)
            res += (exp_lin(resp_tmp[j]) - exp_lin(resp_reduced[j]))*(exp_lin(resp_tmp[j]) - exp_lin(resp_reduced[j]))/(k_s.R[k_s.m->forward_size*j + j]*k_s.R[k_s.m->forward_size*j + j]);

         res = res/(2*simple_model.num_freqs + simple_model.num_channels);
         std::cout << "-------------- " << res << "\n";

         if(res < 1) break;
         if(res > 1.1 * res_old) 
         {
            memcpy(k_s.m->params, old_pars, k_s.m->num_pars*sizeof(double));
            res = res_old;
            break;
         }

         res_old = res;
         //memcpy(k_s.S, S0_s, k_s.m->num_pars * k_s.m->num_pars * sizeof(double));
      }
      std::cout << "REDUCED MODEL FITTED\n";
      simple_model.print_model();

      std::cout << "\nmes     ";
      for(i=0; i<k_s.m->forward_size; i++)
         std::cout << exp_lin(resp_reduced[i]) << "  ";
         
      std::cout << "\nRESP: ";
      for(j=0; j<k_s.m->forward_size; j++)
         std::cout << exp_lin(resp_tmp[j]) << "  "; 

      model.hor_dist = parsed_vals[hd_pos];
      model.ver_dist = parsed_vals[vd_pos];
      model.alt = parsed_vals[alt_pos];

      std::cout << "FITTING FULL MODEL...\n";

      if(started == 1) 
         for(i=0; i<model.num_layers; i++)
            model.rhos[i] = simple_model.rhos[0];
      else
      {
         weight = std::min(1./res_prev, 1.);
         for(i=0; i<model.num_layers; i++)
            model.rhos[i] = weight*model.rhos[i] + (1-weight)*simple_model.rhos[0];
      }

      conv.raw_to_conv();

      k.m->response(resp_tmp);
      res = 0;
      for(j=0; j<k.m->forward_size; j++)
         res += (exp_lin(resp_tmp[j]) - exp_lin(resp[j]))*(exp_lin(resp_tmp[j]) - exp_lin(resp[j]))/(R_0[k.m->forward_size*j + j]*R_0[k.m->forward_size*j + j]);
      res = res/(2*model.num_freqs + model.num_channels);
      res = sqrt(res);

      std::cout << "\nmes     ";
      for(i=0; i<k.m->forward_size; i++)
         std::cout << exp_lin(resp[i]) << "  ";
         
      std::cout << "\nRESP: ";
      for(j=0; j<k.m->forward_size; j++)
         std::cout << exp_lin(resp_tmp[j]) << "  "; 
            
      std::cout << "\nRES: " << res << "\n";
      std::cout << " \n\n\n\n";

      // adjust for log-lin
      for(j=0; j<k.m->forward_size; j++) 
      {
         if (abs(resp[j]) > 5*R_0[j*k.m->forward_size + j])
            k.R[j*k.m->forward_size + j] = R_0[j*k.m->forward_size + j]/resp[j];
         else
            k.R[j*k.m->forward_size + j] = R_0[j*k.m->forward_size + j];
      }

      memcpy(k.S, S0, k.m->num_pars * k.m->num_pars * sizeof(double));
      for(it=0; it<1; it++)
      {
         res_old = res;
         memcpy(best_params, k.m->params, k.m->num_pars*sizeof(double));

         k.compute_sigma();  
         k.get_matrices();  
         k.get_update_vector(upd_vec, upd_cov, resp_reduced);      
         factor = k.m->update(upd_vec, resp_reduced, res, k.R);
         k.update_cov(upd_cov, factor);
         
         k.m->response(resp_tmp);
         
         res = 0;
         for(j=0; j<k.m->forward_size; j++)
            res += (exp_lin(resp_tmp[j]) - exp_lin(resp[j]))*(exp_lin(resp_tmp[j]) - exp_lin(resp[j]))/(R_0[k.m->forward_size*j + j]*R_0[k.m->forward_size*j + j]);

         res = res/(2*model.num_freqs + model.num_channels);
         res = sqrt(res);
         std::cout << "-------------- " << res << "\n";

         if(res<STOP_VAL)
            break;
          
         if(res > 1.1 * res_old) 
         {
            memcpy(k.m->params, best_params, k.m->num_pars*sizeof(double));
            res = res_old;
            break;
         }
         res_old = res;
         //memcpy(k.S, S0, k.m->num_pars * k.m->num_pars * sizeof(double));
      }

      std::cout << "FULL MODEL FITTED\n";
      std::cout << "RES ----------- " << res << "\n";
       std::cout << "\nmes     ";
      for(i=0; i<k.m->forward_size; i++)
         std::cout << exp_lin(resp[i]) << "  ";
         
      std::cout << "\nRESP: ";
      for(j=0; j<k.m->forward_size; j++)
         std::cout << exp_lin(resp_tmp[j]) << "  "; 
      res_prev = res;
      model.print_model();

      output << time << " ";
      output << res << " ";
      model.print_to_file(output);

      output << simple_model.rhos[0] << " ";

      var_mes(S0, k.S, uncert_buf, k.m->num_pars);

      for(i=0; i<k.m->num_pars; i++) 
         output << uncert_buf[i] << " ";

      output << "\n";
   }

   return 0;
}