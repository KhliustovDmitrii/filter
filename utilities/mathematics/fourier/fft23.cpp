#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "fft23.h"

#define SQRT3D2 .8660254037844386

void inv_ind(FFT *fft) 
{
    int i,j,k,l,m,num = 1,offs,order = 0;
    int deo = !(fft->order2);
    for(i=0;i<fft->n;i++)
   fft->inv_index[i] = i;

    for(;order<fft->order3-deo;order++) {
   l = (order&1);
   k = !l;
   if(k) k*=fft->n;
   else  l*=fft->n;
   m = fft->n/num - 1;
   for(i=0;i<num;i++) {
       offs = i*(m+1);
       for(j=0;j<m;j++) {
      fft->inv_index[j+offs+k] = fft->inv_index[(j*3)%m+offs+l];
       }
       fft->inv_index[m+offs+k] = fft->inv_index[m+offs+l];
   }
   num *= 3;
    }

 
    for(;order<fft->order2+fft->order3-1;order++) {
   l = (order&1);
   k = !l;
   if(k) k*=fft->n;
   else  l*=fft->n;
   m = fft->n/num - 1;
   for(i=0;i<num;i++) {
       offs = i*(m+1);
       for(j=0;j<m;j++) {
      fft->inv_index[j+offs+k] = fft->inv_index[(j*2)%m+offs+l];
       }
       fft->inv_index[m+offs+k] = fft->inv_index[m+offs+l];
   }
   num *= 2;
    }
    if(k) 
   memcpy(fft->inv_index,&fft->inv_index[k],k*sizeof(int));
}

void init_fft(FFT *fft, int n)
{
    using namespace std::complex_literals;
    int nn = n,order = 1;
    fft->order2 = fft->order3 = 0;
    fft->n = 1;

    while (nn&&(!(nn&1))) {
   nn >>= 1;
   fft->order2++;
   fft->n <<= 1;
    }

    while (nn&&(!(nn%3))) {
   nn/=3;
   fft->order3++;
   fft->n *= 3;
    }

    fft->inv_index = (int *)realloc(fft->inv_index,2*fft->n*sizeof(int));
    fft->wkn= (std::complex<double> *)realloc(fft->wkn,fft->n*sizeof(std::complex<double>));
    fft->xn = (std::complex<double> *)realloc(fft->xn,fft->n*sizeof(std::complex<double>));
    fft->fn = (std::complex<double> *)realloc(fft->fn,fft->n*sizeof(std::complex<double>));
    
    memset(fft->inv_index,0,2*fft->n*sizeof(int));
    memset(fft->wkn,0,fft->n*sizeof(std::complex<double>));
    memset(fft->xn,0,fft->n*sizeof(std::complex<double>));
    memset(fft->fn,0,fft->n*sizeof(std::complex<double>));

    int i,j;

    inv_ind(fft);

    for(i=0;i<fft->n;i++) {
   fft->wkn[i] = std::exp(-2.*i*1i*M_PI/static_cast<double>(fft->n));
    }
}

int wkn(int k, int n_c, FFT *fft)
{
return (k*(fft->n/n_c));
}

FFT* FFT_new(int w)
{
    FFT *f;
    f=(FFT*)malloc(sizeof(FFT));
    memset(f,0,sizeof(FFT));
    init_fft(f,w);
    return f;
}

void FFT_free(FFT *f)
{
    free(f);
}

void FFT_free_full(FFT *f)
{
    free(f->inv_index);
    free(f->wkn);
    free(f->xn);
    free(f->fn);
    //free(f);
}

void fft_pro(FFT *fft,int inv)
{
    using namespace std::complex_literals;
    int n, nd2, k, mpnd2, mpnd3, m, order;
    std::complex<double> temp,temp1,temp2;

for(m=0;m<fft->n;m++)
   fft->fn[m] = fft->xn[fft->inv_index[m]];

for(n=1,order=0;order<fft->order2;order++) {
    nd2 = n, n+= n;
    for(k=0;k<nd2;k++) {
   int ind = wkn(k,n,fft);
   for (m=k; m<fft->n; m+=n) {
       mpnd2 = m + nd2;
        temp = fft->wkn[ind] * fft->fn[mpnd2];
        fft->fn[mpnd2] = fft->fn[m] - temp;
       fft->fn[m]     = fft->fn[m] + temp;
   }
    }
    
}

for(;order<fft->order2+fft->order3;order++) {
    nd2 = n, n+= 2*n;
    for(k=0;k<nd2;k++) {
   int ind  = wkn(k,n,fft);
   int ind1 = wkn(2*k,n,fft);
   for (m=k; m<fft->n; m+=n) {
       mpnd2 = m + 2*nd2;
       mpnd3 = m + nd2;
       temp  = fft->wkn[ind] * fft->fn[mpnd3];
       temp1 = fft->wkn[ind1]* fft->fn[mpnd2];
       temp2 = temp + temp1;   // sum
       temp  = temp - temp1;  
       temp1 = 1i*(real(temp)*SQRT3D2) -   (imag(temp)*SQRT3D2);  // i*diff*sqrt3/2
       temp  =   (.5*real(temp2))     + 1i*(.5*imag(temp2))    ;  // half sum
       fft->fn[mpnd3] = fft->fn[m] - temp - temp1;
       fft->fn[mpnd2] = fft->fn[m] - temp + temp1;
       fft->fn[m]     = fft->fn[m] + temp2;
   }
    }
}

if(!inv) {
    for(m=1;m<fft->n/2;m++)
        {
        fft->fn[m] += conj(fft->fn[fft->n - m]);
        fft->fn[fft->n - m] = 0;
        fft->fn[m] = conj(fft->fn[m])/static_cast<double>(fft->n);
        }
    fft->fn[0] = 0;
    fft->fn[fft->n/2] = 0;
} else {
    ;//for(m=0;m<fft->n;m++)
     //   fft->fn[m] /= fft->n;
}

}


