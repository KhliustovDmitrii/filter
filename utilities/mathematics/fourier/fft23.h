#ifndef FFT23_H
#define FFT23_H

#include <complex>

namespace filter::math
{
typedef struct
{
    int n;
    int order2;
    int order3;
    int *inv_index;
    std::complex<double> *wkn;
    std::complex<double> *xn;
    std::complex<double> *fn;
} FFT;

FFT *FFT_new(int);
void FFT_free(FFT *);
void FFT_free_full(FFT *);
void fft_pro(FFT *,int);
void init_fft(FFT *fft, int n);
void fft_set_zero(FFT *);

}; // filter::math
#endif
