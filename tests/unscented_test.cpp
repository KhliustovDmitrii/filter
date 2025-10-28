// g++ unscented_test.cpp ../filters/kalman_unscented.cpp ../models/aem_model_raw.cpp ../converters/aem_converter.cpp ../fft23.cpp ../text_utils.cpp ../math_utils.cpp -o forward

#include <iostream>
#include <fstream>

#include "../models/aem_model_raw.h"
#include "../converters/aem_converter.h"
#include "../filters/kalman_unscented.h"
#include "../fft23.h"
#include "../text_utils.h"
#include "../math_utils.h"


#define PROCESS_TEST 1
#define STOP_VAL 1.


int main(int argc, char *argv[])
{
   int i, j, it, l;
   double res, res_old;

   std::ifstream data_file(argv[1]);
   std::ifstream wf("waveform.XYZ");
   std::ofstream out(argv[2]);
   std::ofstream out_var(argv[3]);

   std::string buf;

   int num_freqs = 10;
   double *freqs_in = new double[num_freqs];

   double chans_in[15] = {1, 2, 3, 5, 7, 12, 19, 30, 49, 79, 128, 207, 336, 545, 884};
   int num_pol_layers = 1;
   int *pol_lays = new int[num_pol_layers+1];
   double bfr = 77.16;
   double *waveform = new double[2592];

   int num_layers = 1;

   pol_lays[0] = 0;

   for(i=0; i<num_freqs; i++) freqs_in[i] = (i+1)*100;

   aem_model_raw model = aem_model_raw(num_layers, num_pol_layers, num_freqs, 14, 100, 2592, 199, freqs_in, chans_in, pol_lays);

   for(i=1;i<2*model.num_freqs_fulltime;i+=2)
      model.freqs_fd_fulltime[(i-1)/2] = bfr*i;

    for(i=0; i<num_layers-1; i++)
        model.depths[i] = (i+1)*10;

    for(i=0; i<num_layers; i++)
        model.rhos[i] = (i+1)*90;

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

   int diff_ac = 0;
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

   kalman_unscented k = kalman_unscented(&conv);

   k.R = new double[k.m->forward_size*k.m->forward_size];
   k.S = new double[k.m->num_pars*k.m->num_pars];

   for(i=0; i<k.m->num_pars; i++)
   {
      for(j=0; j<k.m->num_pars; j++)
      {
         if(i==j)
            k.S[k.m->num_pars*i+j] = 0.1;
         else
            k.S[k.m->num_pars*i+j] = 0;
      }
   }

   for(i=0; i<k.m->forward_size; i++)
   {
      for(j=0; j<k.m->forward_size; j++)
      {
         if(i==j)
            k.R[k.m->forward_size*i+j] = 0.1;
         else
            k.R[k.m->forward_size*i+j] = 0;
      }
   }

   for(j=2*model.num_freqs; j<2*model.num_freqs+model.num_channels; j++)
      k.R[j*k.m->forward_size+j] = 50000;

   double resp[k.m->forward_size];
   double resp_tmp[k.m->forward_size];
   double resp_parsed[k.m->forward_size];
   double best_params[k.m->num_pars];
   double upd_vec[k.m->num_pars];
   double upd_cov[k.m->num_pars*k.m->num_pars];
   double factor;

   while(getline(data_file, buf))
   {
      parse(buf, resp_parsed);

      for(j=0; j<model.num_freqs; j++) resp[2*j] = resp_parsed[j + model.num_freqs];
      for(j=0; j<model.num_freqs-1; j++) resp[2*j+1] = resp_parsed[j+1]-resp_parsed[j];
      resp[2*model.num_freqs - 1] = resp[2*model.num_freqs - 3];
      for(j=2*model.num_freqs; j<2*model.num_freqs+model.num_channels; j++) resp[j] = resp_parsed[j];
      for(j=0; j<k.m->forward_size; j++) resp[j] = log_lin(resp[j]);

      k.m->response(resp_tmp);
         
      res = 0;
      for(j=0; j<k.m->forward_size; j++)
         res += (exp_lin(resp_tmp[j]) - exp_lin(resp[j]))*(exp_lin(resp_tmp[j]) - exp_lin(resp[j]))/(k.R[k.m->forward_size*j + j]*k.R[k.m->forward_size*j + j]);
      res = res/(2*model.num_freqs + model.num_channels);
   

      std::cout << "\nmes     ";
      for(i=0; i<k.m->forward_size; i++)
         std::cout << resp[i] << "  ";
         
      std::cout << "\nRESP: ";
      for(j=0; j<k.m->forward_size; j++)
         std::cout << resp_tmp[j] << "  "; 
            
      std::cout << "\nRES: " << res << "\n";
      std::cout << "\nPAR " << "\n";
         for(j=0; j<k.m->num_pars; j++)
            std::cout << k.m->params[j] << "  ";
      std::cout << " \n\n\n\n";


      for(it=0; it<100; it++)
      {
         res_old = res;
         memcpy(best_params, k.m->params, k.m->num_pars*sizeof(double));
         
         k.compute_sigma();  
         k.get_matrices();  
         k.get_update_vector(upd_vec, upd_cov, resp);    
         for(i=0; i<k.m->num_pars; i++) upd_vec[i] = -upd_vec[i];  
         factor = k.m->update(upd_vec, resp, res, k.R);
         k.update_cov(upd_cov, factor);
         
         k.m->response(resp_tmp);
         
         res = 0;
         for(j=0; j<k.m->forward_size; j++)
            res += (exp_lin(resp_tmp[j]) - exp_lin(resp[j]))*(exp_lin(resp_tmp[j]) - exp_lin(resp[j]))/(k.R[k.m->forward_size*j + j]*k.R[k.m->forward_size*j + j]);

         res = res/(2*model.num_freqs + model.num_channels);
         /*
         std::cout << "\nSTEP " << it << "  ";
         for(j=0; j<k.m->num_pars; j++)
            std::cout << k.m->params[j] << "  ";
            
               
         std::cout << "\nRESP ";
         for(j=0; j<k.m->forward_size; j++)
            std::cout << resp_tmp[j] << "  ";   
            
         */     
         std::cout << "\nRES " << res << " \n";
         

         //out << res << " ";
         
         /*
         std::cout << "\nUPD ";
         for(j=0; j<k.m->num_pars; j++)
            std::cout << upd_vec[j] << "  ";
         */    
         std::cout << " \n\n";
         
         if(res<STOP_VAL)
            break;
          

         /*
         if(res>1.05*res_old && it > 1)
         {
            break;
         }
         */
         /*
         for(j=0; j<k.m->num_pars; j++)
         {
            for(l=0; l<k.m->num_pars; l++)
            {
               if(j==l)
                  k.S[k.m->num_pars*j+l] = 1;
               else
                  k.S[k.m->num_pars*j+l] = 0;
            }
         }
         */
         model.print_model();
      }
      memcpy(k.m->params, best_params, k.m->num_pars*sizeof(double));
      model.print_model();
   }

   for(j=0; j<k.m->num_pars; j++)
   {
      for(l=0; l<k.m->num_pars; l++)
      {
         out_var << k.S[k.m->num_pars*j+l] << " ";
      }
      out_var << "\n";
   }

   return 0;
}
