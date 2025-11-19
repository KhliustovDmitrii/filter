#ifndef _EQUATOR_H_
#define _EQUATOR_H_
// main model for airborne electromagnetics
// transmitter and receiver are coils parallel to the ground

#include <vector>
#include <complex>

#include "../../utilities/mathematics/fourier/fft23.h"

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
            int nff = 100, int sl = 2592, int bc = 199, double bfr = 77.16,
            std::vector<std::complex<double>> spec) :
            num_layers(nl), num_pol_layers(pi.size()), pol_inds(pi), 
            num_freqs(fr.size()), num_channels(ch.size() - 1),
            depths(std::vector<double>(0, nl-1)), 
            rhos(std::vector<double>(rho, nl)),
            cole_rho(std::vector<double>(rho, pi.size())),
            cole_tau(std::vector<double>(0.001, pi.size())),
            cole_c(std::vector<double>(0.5, pi.size())),
            altitude_correction(std::vector<double>(0, 1)), 
            freqs(fr), chans(ch),
            freq_re(std::vector<double>(0, fr.size())),
            freq_im(std::vector<double>(0, fr.size())),
            td_chan(std::vector<double>(0, ch.size()-1)),
            hor_dist(0), ver_dist(0), alt(0), prim_field(pf),
            num_freqs_fulltime(nff), spec_len(sl), base_chan(bc),
            freqs_fd_fulltime(std::vector<double>(0, nff)),
            impulse_spec(spec)
    {
        int i;
        double thick_0 = 4;

        for(i=0; i<num_layers; i++)
        {
            depths[i] = thick_0;
            thick_0 = thick_0 * 1.1085;
        }
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