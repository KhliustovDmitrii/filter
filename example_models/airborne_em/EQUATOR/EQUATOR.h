#ifndef _EQUATOR_H_
#define _EQUATOR_H_
// main model for airborne electromagnetics
// transmitter and receiver are coils parallel to the ground

class EQUATOR
{
public:

    // dimensionality
    int num_layers;
    int num_pol_layers;
    int pol_inds[num_pol_layers];

    int num_freqs;
    int num_channels;
    int num_freqs_fulltime;
    int spec_len;
    int base_chan;

    // parameters
    double depths[num_layers-1];
    double rhos[num_layers];
    double cole[3*num_pol_layers];

    double freqs[num_freqs];
    double chans[num_channels+1];

    double altitude_correction;

    double freqs_fd_fulltime[num_freqs_fulltime];

    // response
    double freq_re[num_freqs];
    double freq_im[num_freqs];
    double td_chan[num_channels];

    // geometry
    double hor_dist;
    double ver_dist;
    double alt;
    double prim_field;

    // auxillary

    // frequency dependent resistivities
    std::complex<double> rhos_fd[num_layers];

    // for response calculation
    std::complex<double> resp_fd[num_freqs];
    std::complex<double> resp_fd_fulltime[num_freqs_fulltime];


    std::complex<double> impulse_spec[spec_len];
    FFT fft;

    // forward problem
    void forward();

    // pretty printing
    void print_model();
    void print_to_file(std::ofstream &out);

    // auxillary functions for computations
    std::complex<double> PartSum(double n0, double hh, double r, std::complex<double> n1, std::complex<double> Imp);
    std::complex<double> Impedance(double n0, double om);
    std::complex<double> integral(double hh, double r, double f);
    std::complex<double> ImHz(double r,double z,double f);
    void fd_forward(std::complex<double> *dest);
    void td_forward(double *dest);
}

#endif