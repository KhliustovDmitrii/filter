#include "core/filter/kalman_extended/kalman_extended.h"
#include "core/filter/kalman_unscented/kalman_unscented.h"
#include "core/input_filter/averaging_filter/averaging_filter.h"
#include "core/regularizer/boundary_projector/boundary_projector.h"
#include "core/updater/decay_updater/decay_updater.h"

#include "model/EQUATOR/EQUATOR.h"
#include "model/EQUATOR_C/EQUATOR_C.h"
#include "data_loader/EQUATOR_data_loader.h"

#include "utilities/parsing/simple_parser/simple_parser.h"
#include "utilities/mathematics/special_functions/special_functions.h"

#include "types/workspace/filter_workspace.h"

#include <iostream>
#include <vector>

// Example of using the filter library for processing airborne electromagnetics data
int main(int argc, char *argv[])
{
  // ------------ OUTPUT FILE
  std::ofstream out(argv[2]);

  // ------------ WAVEFORM
  std::ifstream wf("waveform.XYZ");
  std::vector<std::complex<double>> waveform;
  int i;
  

  for(i=0; i<2592; i++)
  {
    double temp;
    wf >> temp;
    waveform.push_back(temp);
  }

  filter::math::Fourier wf_fft(2592);
  
  wf_fft.set_xn(waveform);
  wf_fft.compute(0);

  // ------------ MODEL

  // full
  int num_layers = 25;
  int num_freqs = 15;
  int num_channels = 14;
  std::vector<int> pol_inds = {};
  double rho_ini = 100;

  std::vector<double> freqs = {77.1, 231.5, 385.8, 540.1, 694.4, 1003.1, 1157.41, 1311.7, 1466, 1620.4, 1774.7,
             1929, 2083.3, 2237.65, 2854.94, 3009.26, 3163.58, 3317.9, 3472.22, 6172.84, 13503.09};

  //std::vector<double> chans = {1, 2, 3, 5, 7, 12, 19, 30, 49, 79, 128, 207, 336, 545, 884};
  std::vector<double> chans = {0};

  double prim_field = filter::math::primField(32, 25)/10000;

  std::vector<std::complex<double>> spec = wf_fft.get_fn();

  filter::examples::EQUATOR full_model(num_layers, pol_inds, rho_ini, freqs, chans, prim_field, spec);

  // reduced
  filter::examples::EQUATOR reduced_model(1, {}, 100, {77.1, 231.5}, {0}, prim_field, {1});

  // ------------ ADAPTER

  // full
  std::vector<double> noise_fd = {0.2, 0.12, 0.09, 0.13, 0.39, 0.5, 0.58, 0.52, 0.36, 0.46, 0.15, 0.31, 1.01, 0.56, 0.63,
             0.76, 0.44, 0.23, 0.44, 0.52, 1.11};
  //std::vector<double> noise_td = {71000, 71000, 70000, 66000, 53000, 33000, 33000, 23000, 19000, 30000, 18000, 16000, 10000, 7000};
  std::vector<double> noise_td = {};
   
  std::vector<double> weights = {};

  for(i=0; i<noise_fd.size(); i++) weights.push_back(noise_fd[i]);
  for(i=0; i<noise_fd.size(); i++) weights.push_back(noise_fd[i]);
  for(i=0; i<noise_td.size(); i++) weights.push_back(noise_td[i]);

  filter::examples::EQUATOR_C adapter(full_model, weights);

  for(i=0; i<full_model.num_layers; i++)
    adapter.add_param<filter::math::ExpLin, filter::math::LogLin>(full_model.rhos, i);

  adapter.add_param<filter::math::Lin, filter::math::Lin>(full_model.altitude_correction, 0);

  // reduced
  filter::examples::EQUATOR_C adapter_red(reduced_model, {0.2, 1000, 0.2, 0.2});


  for(i=0; i<reduced_model.num_layers; i++)
    adapter_red.add_param<filter::math::ExpLin, filter::math::LogLin>(reduced_model.rhos, i);

  // ------------ READER
  char sep = ' ';
  char dec = '.';
  int time_size = 16;
  std::vector<std::string> drops;
  drops.push_back("/");
  drops.push_back("L");

  filter::io::Simple_Parser parser(std::string(argv[1]), sep, dec, time_size, drops);

  // ------------ INPUT FILTER
  filter::components::Averaging_Filter input_filter(4, 256);

  // ------------ DATA LOADER

  // full
  filter::examples::EQUATOR_data_loader data_loader(adapter, 0, 1, 2, 3, 24, 45, -1, 1, 1);

  // reduced
  filter::examples::EQUATOR_data_loader data_loader_red(adapter_red, 0, 1, 2, 3, 24, 45, -1, 1, 1);

  // ------------ FILTER - EXTENDED

  // full
  filter::Kalman_Unscented filter_ext(adapter);
  std::vector<double> R = filter_ext.get_R();
  std::vector<double> S = filter_ext.get_S();

  double err_ini = 0.3;
  double cor_ini = 0.1;
  i = adapter.num_pars-1;
  S[i+adapter.num_pars*i] = 0.3; // altitude correction

  for(i=adapter.num_pars-2;i>=0;i--) 
  {
    if(i==adapter.num_pars-2)//Bottom layer resisitivity
     S[i+adapter.num_pars*i] = err_ini;
    else //Other resisitivities
    {
     S[i+1+adapter.num_pars*i] = cor_ini*err_ini/S[i+1+adapter.num_pars*(i+1)];
     S[i+adapter.num_pars*i] = sqrt(err_ini*err_ini-S[i+1+adapter.num_pars*i]*S[i+1+adapter.num_pars*i]);
    }
  }

  // real component
  for(i=0; i<full_model.num_freqs-1; i++)
    R[adapter.forward_size*i+i] = sqrt(noise_fd[i]*noise_fd[i] + noise_fd[i+1]*noise_fd[i+1]);

  // higher freq excluded
  R[adapter.forward_size*i + i] = 1000;
  i++;

  // imaginary component
  int pos = i;
  for(i=pos; i<full_model.num_freqs+pos; i++)
    R[adapter.forward_size*i+i] = noise_fd[i-pos];

  // td
  pos = i;
  for(i=pos; i<pos+full_model.num_channels; i++)
    R[adapter.forward_size*i+i] = noise_td[i-pos];

  filter_ext.set_R(R);
  filter_ext.set_S(S);

  auto extended_ws_full = filter_ext.allocate_workspace();

  // reduced
  filter::Kalman_Unscented filter_ext_red(adapter_red);
  R = filter_ext_red.get_R();
  S = filter_ext_red.get_S();

  S[0] = 0.3;
  R[0] = 0.2;
  R[5] = 1000;
  R[10] = 0.2;
  R[15] = 1000;

  filter_ext_red.set_R(R);
  filter_ext_red.set_S(S);

  auto extended_ws_red = filter_ext_red.allocate_workspace();
  // ------------ UPDATER
  filter::Decay_Updater updater(adapter, 3, 0.5);
  filter::Decay_Updater updater_red(adapter_red, 3, 0.5);

  // ------------ REGULARIZER
  // TODO

  int ret_code;
  std::vector<double> read_result(256, 0);
  std::vector<double> avg_result(256, 0);
  std::string time;
  std::vector<double> measurements(256, 0);
  std::vector<double> measurements_red(256, 0);
  std::vector<double> upd_vec(adapter.num_pars, 0);
  std::vector<double> upd_cov(adapter.num_pars*adapter.num_pars, 0);

  std::vector<double> response(adapter.forward_size, 0);

  double residual, residual_min;

  std::vector<double> best_params(adapter.num_pars, 0);
  double best_full_residual = 1000;
  
  while((ret_code = parser.read(read_result, time)) != -1)
  {
    // skip the line
    if(ret_code == 1) continue;

    ret_code = input_filter.add_data(read_result, avg_result);

    // more data to add
    if(ret_code != 0) continue;

    // fit reduced model
    data_loader_red.load(avg_result);
    data_loader_red.get_measurements(measurements_red);

    adapter_red.response(response);
    residual_min = adapter_red.residual(measurements_red, response);

    for(i=0; i<10; i++)
    {
      filter_ext_red.get_update(measurements_red, upd_vec, upd_cov, *extended_ws_red);
      updater_red.update(upd_vec, measurements_red);

      adapter_red.response(response);
      residual = adapter_red.residual(measurements_red, response);
      std::cout << "------   " << residual << std::endl;

      if(residual < 1 || residual > residual_min) break;

      residual_min = residual;
    }

    reduced_model.print_model();

    // fit full model
    data_loader.load(avg_result);
    data_loader.get_measurements(measurements);

    // weight full and reduced model
    double weight = std::min(1., 1./best_full_residual);
    for(i=0; i<full_model.num_layers; i++)
    {
      full_model.rhos[i] = weight*full_model.rhos[i] + (1-weight)*reduced_model.rhos[0];
    }

    adapter.response(response);
    residual_min = adapter.residual(measurements, response);

    for(i=0; i<15; i++)
    {
      filter_ext.get_update(measurements, upd_vec, upd_cov, *extended_ws_full);
      updater.update(upd_vec, measurements);

      adapter.response(response);
      residual = adapter.residual(measurements, response);
      std::cout << "++++++   " << residual << std::endl;

      residual_min = std::min(residual, residual_min);

      if(residual < 1 || residual > residual_min) break;
      else // save parameters
      {
        residual_min = residual;
        for(int par_num = 0; par_num < adapter.num_pars; par_num++)
        {
          best_params[i] = adapter.get_param(par_num);
        }
      }
    }

    if(residual > residual_min) // restore the best model
    {
      for(int par_num = 0; par_num < adapter.num_pars; par_num++)
      {
        adapter.set_param(par_num, best_params[i]);
      }
    }

    best_full_residual = residual_min;

    full_model.print_model();

    out << time << " ";
    for(int rho_num = 0; rho_num < full_model.num_layers; rho_num++)
    {
      out << full_model.rhos[rho_num] << " ";
    }

    out << full_model.altitude_correction[0] << " ";

    double cum_depth = 0;
    for(int d_num = 0; d_num < full_model.num_layers-1; d_num++)
    {
      out << cum_depth + full_model.depths[d_num]/2 << " ";
      cum_depth += full_model.depths[d_num];
    }

    out << cum_depth + full_model.depths[full_model.num_layers - 2] << " ";
    out << reduced_model.rhos[0] << std::endl;
  }

  return 0;
}