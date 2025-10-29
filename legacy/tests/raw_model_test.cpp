//g++ raw_model_test.cpp ../filters/kalman_unscented.cpp ../models/aem_model_raw.cpp ../converters/aem_converter.cpp ../fft23.cpp ../text_utils.cpp ../math_utils.cpp -o filter

#include <iostream>
#include <fstream>
#include "../models/aem_model_raw.h"
#include "../fft23.h"

int main(int argc, char *argv[])
{
    std::ofstream out(argv[1]);
    std::ifstream wf("waveform.XYZ");

    int i;
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

    model.forward();

    for(i=0; i<num_freqs; i++) out << model.freq_re[i] << " ";
    for(i=0; i<num_freqs; i++) out << model.freq_im[i] << " ";
    for(i=0; i<14; i++) out << std::fixed << model.td_chan[i] << " ";
    out << "\n";

    return 0;
}

