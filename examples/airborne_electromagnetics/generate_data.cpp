#include "core/filter/kalman_extended/kalman_extended.h"
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

// Generate synthetic data
int main(int argc, char *argv[])
{
  // ------------ WAVEFORM
  std::ifstream wf("waveform.XYZ");
  std::vector<double> waveform(2592, 0);
  int i;
  

  for(i=0; i<2592; i++)
  {
    double temp;
    wf >> temp;
    waveform.push_back(temp);
  }

  filter::math::FFT imp_fft;
  memset(&imp_fft, 0, sizeof(imp_fft));
  
  filter::math::init_fft(&imp_fft,2592);
  for(i=0;i<2592;i++)
    imp_fft.xn[i] = waveform[i];
    
  filter::math::fft_pro(&imp_fft,0);

  // ------------ MODEL

  // full
  int num_layers = 5;
  int num_freqs = 15;
  int num_channels = 14;
  std::vector<int> pol_inds = {};
  double rho_ini = 10;

  std::vector<double> freqs = {77.1, 231.5, 385.8, 540.1, 694.4, 1003.1, 1157.41, 1311.7, 1466, 1620.4, 1774.7,
             1929, 2083.3, 2237.65, 2854.94, 3009.26, 3163.58, 3317.9, 3472.22, 6172.84, 13503.09};

  //std::vector<double> chans = {1, 2, 3, 5, 7, 12, 19, 30, 49, 79, 128, 207, 336, 545, 884};
  std::vector<double> chans = {0};

  double prim_field = filter::math::primField(32, 25)/10000;

  std::vector<std::complex<double>> spec;
  for(i=0;i<2592;i++)
    spec.push_back(imp_fft.fn[i]);

  filter::examples::EQUATOR full_model(num_layers, pol_inds, rho_ini, freqs, chans, prim_field, spec);

  // GENERATE DATA

  std::ofstream out("synthetic_data.XYZ");

  for(int num_line = 0; num_line<10; num_line++)
  {
    out << num_line << "     ";

    full_model.hor_dist = 25;
    full_model.ver_dist = 25;
    full_model.alt = 25;
    full_model.forward();

    out << full_model.hor_dist << "  " << full_model.ver_dist << "  " << full_model.alt << "  ";

    for(i=0; i<full_model.num_freqs; i++)
        out << full_model.freq_re[i] << "  ";

    for(i=0; i<full_model.num_freqs; i++)
        out << full_model.freq_im[i] << "  ";

    for(i=0; i<full_model.num_channels; i++)
        out << full_model.td_chan[i] << "  ";

    out << std::endl;
  }
  
  return 0;
}