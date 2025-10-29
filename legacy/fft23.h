#ifndef _FFT23_H_

#define _FFT23_H_

#include <complex>


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


#endif
