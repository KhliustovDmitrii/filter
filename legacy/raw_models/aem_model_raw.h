#ifndef _AEM_MODEL_RAW_H_
#define _AEM_MODEL_RAW_H_

#include <complex>
#include <cstring>
#include <fstream>
#include "model.h"
#include "../fft23.h"

#define mu0 (M_PI*.4e-6)

class aem_model_raw
{
public:

   // dimensionality
   int num_layers;
   int num_pol_layers;
   int *pol_inds;

   int num_freqs;
   int num_channels;
   int num_freqs_fulltime;
   int spec_len;
   int base_chan;

   // parameters
   double *depths;
   double *rhos;
   double *cole;

   double *freqs;
   double *chans;

   double altitude_correction;

   double *freqs_fd_fulltime;

   // response
   double *freq_re;
   double *freq_im;
   double *td_chan;

   // geometry
   double hor_dist;
   double ver_dist;
   double alt;
   double prim_field;

   // auxillary

   // frequency dependent resistivities
   std::complex<double> *rhos_fd;

   // for response calculation
   std::complex<double> *resp_fd;
   std::complex<double> *resp_fd_fulltime;


   std::complex<double> *impulse_spec;
   FFT fft;

   // only init variables that do not change during execution
   aem_model_raw(int nl, int npl, int nf, int nc, int nff, int sl, int bchn, 
                                                 double *freqs_in, 
                                                 double *chans_in,
                                                 int *pol_lays)
   {
      // input dimension
      num_layers = nl;
      num_pol_layers = npl;
      num_freqs_fulltime = nff;
      spec_len = sl;
      base_chan = bchn;

      pol_inds = new int[num_pol_layers];
      std::memcpy(pol_inds, pol_lays, num_pol_layers*sizeof(int));

      // parameters
      depths = new double[num_layers-1];
      rhos = new double[num_layers];
      cole = new double[3*num_pol_layers];

      // parameter dimension
      num_freqs = nf;
      num_channels = nc;

      freqs = new double[num_freqs];
      chans = new double[num_channels+1];
      std::memcpy(freqs, freqs_in, num_freqs*sizeof(double));
      std::memcpy(chans, chans_in, (num_channels+1)*sizeof(double));

      // output dimension
      freq_re = new double[num_freqs];
      freq_im = new double[num_freqs];
      td_chan = new double[num_channels];

      // auxillary
      rhos_fd = new std::complex<double>[num_layers];
      resp_fd = new std::complex<double>[num_freqs];
      resp_fd_fulltime = new std::complex<double>[num_freqs_fulltime];
      impulse_spec = new std::complex<double>[spec_len];

      freqs_fd_fulltime = new double[num_freqs_fulltime];
   }

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
   void fd_forward(std::complex<double> *dest, double *freqs, int num_freqs);
   void td_forward(double *dest);
};

#endif
