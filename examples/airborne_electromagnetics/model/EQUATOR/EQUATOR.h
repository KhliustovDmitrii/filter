#ifndef _EQUATOR_H_
#define _EQUATOR_H_
// main model for airborne electromagnetics
// transmitter and receiver are coils parallel to the ground

#include <vector>
#include <complex>

#include "utilities/mathematics/fourier/fft23.h"

const double mu0 = M_PI*.4e-6;

class EQUATOR
{
public:

    // dimensionality
    int num_layers;
    int num_pol_layers;
    std::vector<int> pol_inds;

    int num_freqs;
    int num_channels;

    // parameters

    // differentiable
    std::vector<double> depths, rhos;
    std::vector<double> cole_rho, cole_tau, cole_c;
    std::vector<double> altitude_correction;
    
    // fixed
    std::vector<double> freqs, chans;

    // response
    std::vector<double> freq_re, freq_im, td_chan;

    // geometry
    double hor_dist;
    double ver_dist;
    double alt;
    double prim_field;

    // forward problem
    void forward();

    // pretty printing
    void print_model();
    void print_to_file(std::ofstream &out);

    // init the model to most often used form - n layer half space of uniform resistivity
    // with lay_thick[i] = 1.1085 * lay_thick[i-1]
    EQUATOR(int nl, std::vector<int> pi, double rho, 
            std::vector<double> fr, std::vector<double> ch, double pf,
            std::vector<std::complex<double>> spec,
            int nff = 100, int sl = 2592, int bc = 199, double bfr = 77.16) :
            num_layers(nl), num_pol_layers(pi.size()), pol_inds(pi), 
            num_freqs(fr.size()), num_channels(ch.size() - 1),
            depths(std::vector<double>(nl-1, 0)), 
            rhos(std::vector<double>(nl, rho)),
            cole_rho(std::vector<double>(pi.size(), rho)),
            cole_tau(std::vector<double>(pi.size(), 0.001)),
            cole_c(std::vector<double>(pi.size(), 0.5)),
            altitude_correction(std::vector<double>(1, 0)), 
            freqs(fr), chans(ch),
            freq_re(std::vector<double>(fr.size(), 0)),
            freq_im(std::vector<double>(fr.size(), 0)),
            td_chan(std::vector<double>(ch.size()-1, 0)),
            hor_dist(0), ver_dist(0), alt(0), prim_field(pf),
            num_freqs_fulltime(nff), spec_len(sl), base_chan(bc),
            freqs_fd_fulltime(std::vector<double>(nff, 0)),
            impulse_spec(spec), 
            rhos_fd(std::vector<std::complex<double>>(nl, 0)),
            resp_fd(std::vector<std::complex<double>>(2*fr.size(), 0)),
            resp_fd_fulltime(std::vector<std::complex<double>>(2*nff, 0))
    {
        int i;
        double thick_0 = 4;

        for(i=0; i<num_layers-1; i++)
        {
            depths[i] = thick_0;
            thick_0 = thick_0 * 1.1085;
        }

        for(i=1;i<2*num_freqs_fulltime;i+=2)
            freqs_fd_fulltime[(i-1)/2] = bfr*i;
    };



private:

    // intrinsic model parameters
    int num_freqs_fulltime;
    int spec_len;
    int base_chan;
    std::vector<double> freqs_fd_fulltime; // freqs for full-time td computation

    // auxillary

    // frequency dependent resistivities
    std::vector<std::complex<double>> rhos_fd;

    // for response calculation
    std::vector<std::complex<double>> resp_fd, resp_fd_fulltime;

    std::vector<std::complex<double>> impulse_spec;
    FFT fft;

    // auxillary functions for computations
    std::complex<double> PartSum(double n0, double hh, double r, std::complex<double> n1, std::complex<double> Imp);
    std::complex<double> Impedance(double n0, double om);
    std::complex<double> integral(double hh, double r, double f);
    std::complex<double> ImHz(double r,double z,double f);
    void fd_forward(std::vector<std::complex<double>> &dest, std::vector<double> &frs);
    void td_forward(std::vector<double> &dest);
};

#endif