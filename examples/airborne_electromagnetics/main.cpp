#include "core/filter/kalman_extended/kalman_extended.h"
#include "core/input_filter/averaging_filter/averaging_filter.h"
#include "core/regularizer/boundary_projector/boundary_projector.h"
#include "core/updater/decay_updater/decay_updater.h"

#include "model/EQUATOR/EQUATOR.h"
#include "model/EQUATOR_C/EQUATOR_C.h"
#include "data_loader/EQUATOR_data_loader.h"

#include "utilities/parsing/simple_parser/simple_parser.h"
#include "utilities/mathematics/special_functions/special_functions.h"


// Example of using the filter library for processing airborne electromagnetics data
int main(int argc, char *argv[])
{
    // ------------ WAVEFORM
    std::ifstream wf("waveform.XYZ");
    double *waveform = new double[2592];
    int i;
   
    for(i=0; i<2592; i++) wf >> waveform[i];

    FFT imp_fft;
    memset(&imp_fft, 0, sizeof(imp_fft));
    
    init_fft(&imp_fft,2592);
    for(i=0;i<2592;i++)
        imp_fft.xn[i] = waveform[i];
        
    fft_pro(&imp_fft,0);

    // ------------ MODEL
    int num_layers = 25;
    int num_freqs = 15;
    int num_channels = 14;
    std::vector<int> pol_inds = {};
    double rho_ini = 100;

    std::vector<double> freqs = {77.1, 231.5, 385.8, 540.1, 694.4, 1003.1, 1157.41, 1311.7, 1466, 1620.4, 1774.7,
                       1929, 2083.3, 2237.65, 2854.94, 3009.26, 3163.58, 3317.9, 3472.22, 6172.84, 13503.09};

    std::vector<double> chans = {1, 2, 3, 5, 7, 12, 19, 30, 49, 79, 128, 207, 336, 545, 884};

    double prim_field = primField(32, 25)/10000;

    std::vector<std::complex<double>> spec;
    for(i=0;i<2592;i++)
        spec.push_back(imp_fft.fn[i]);

    EQUATOR full_model(num_layers, pol_inds, rho_ini, freqs, chans, prim_field, spec);

    // ------------ ADAPTER

    std::vector<int> dr, dac, sr, sac;
    std::vector<int> empty = {};
    for(i=0; i<num_layers; i++)
    {
        dr.push_back(i);
        sr.push_back(2);
    }
    dac.push_back(1);
    sac.push_back(0);

    std::vector<double> noise_fd = {0.2, 0.12, 0.09, 0.13, 0.39, 0.5, 0.58, 0.52, 0.36, 0.46, 0.15, 0.31, 1.01, 0.56, 0.63,
                       0.76, 0.44, 0.23, 0.44, 0.52, 1.11};
    std::vector<double> noise_td = {71000, 71000, 70000, 66000, 53000, 33000, 33000, 23000, 19000, 30000, 18000, 16000, 10000, 7000};
   
    std::vector<double> weights = {};

    for(i=0; i<noise_fd.size(); i++) weights.push_back(noise_fd[i]);
    for(i=0; i<noise_fd.size(); i++) weights.push_back(noise_fd[i]);
    for(i=0; i<noise_td.size(); i++) weights.push_back(noise_td[i]);

    EQUATOR_C adapter(&full_model, empty, dr, dac, empty, empty, empty, empty, sr, sac, empty, empty, empty, weights);

    // ------------ READER
    char sep = ' ';
    char dec = '.';
    int time_size = 15;
    std::vector<std::string> drops;
    drops.push_back("/");
    drops.push_back("L");

    Simple_Parser parser(std::string(argv[1]), sep, dec, time_size, drops);

    // ------------ INPUT FILTER
    Averaging_Filter input_filter(4, 256);

    // ------------ DATA LOADER
    EQUATOR_data_loader data_loader(&adapter, 0, 1, 2, 3, 24, 45);

    // ------------ FILTER - EXTENDED
    Kalman_Extended filter_ext(&adapter);
    std::vector<double> R = filter_ext.get_R();
    std::vector<double> S = filter_ext.get_S();

    double err_ini = 0.3;
    double cor_ini = 0.1;
    for(i=adapter.num_pars-1;i>=0;i--) 
    {
      if(i==adapter.num_pars-1)//Bottom layer resisitivity
         S[i+adapter.num_pars*i] = err_ini;
      else //Other resisitivities
      {
         S[i+1+adapter.num_pars*i] = cor_ini*err_ini/S[i+1+adapter.num_pars*(i+1)];
         S[i+adapter.num_pars*i] = sqrt(err_ini*err_ini-S[i+1+adapter.num_pars*i]*S[i+1+adapter.num_pars*i]);
      }
    }

    for(i=0; i<full_model.num_freqs-1; i++)
      R[adapter.forward_size*i+i] = sqrt(noise_fd[i]*noise_fd[i] + noise_fd[i+1]*noise_fd[i+1]);

    R[adapter.forward_size*i + i] = 1000;
    i++;

    int pos = i;
    for(i=pos; i<full_model.num_freqs+pos; i++)
      R[adapter.forward_size*i+i] = noise_fd[i-pos];

    pos = i;
    for(i=pos; i<pos+full_model.num_channels; i++)
      R[adapter.forward_size*i+i] = noise_td[i-pos];

    filter_ext.set_R(R);
    filter_ext.set_S(S);

    // ------------ UPDATER
    Decay_Updater updater(&adapter, 3, 0.5);

    // ------------ REGULARIZER
    // TODO

    int ret_code;
    std::vector<double> read_result(256, 0);
    std::vector<double> avg_result(256, 0);
    std::string time;
    std::vector<double> measurements(256, 0);
    std::vector<double> upd_vec(adapter.num_pars, 0);
    std::vector<double> upd_cov(adapter.num_pars*adapter.num_pars, 0);
    
    while((ret_code = parser.read(read_result, time)) != -1)
    {
        // skip the line
        if(ret_code == 1) continue;

        ret_code = input_filter.add_data(read_result, avg_result);

        // more data to add
        if(ret_code != 0) continue;

        data_loader.load(avg_result);

        data_loader.get_measurements(measurements);

        filter_ext.get_update(measurements, upd_vec, upd_cov);

        updater.update(upd_vec, measurements);

        full_model.print_model();
    }
    return 0;
}