#include <iostream>
#include "aem_model_raw.h"
#include "../math_utils.h"

void aem_model_raw::forward()
{
    int i;

    // --------------------  frequency domain response
    fd_forward(resp_fd, freqs, num_freqs);

    for(i=0; i<num_freqs; i++)
    {
        freq_re[i] = real(resp_fd[i]);
        freq_im[i] = imag(resp_fd[i]);
    }


    // --------------------  time domain response
    if(num_channels > 0)
        td_forward(td_chan);
}

void aem_model_raw::td_forward(double *dest)
{
    using namespace std::complex_literals;
    spline sr, si;
    int i, j;
   
    // fulltime frequency response
    fd_forward(resp_fd_fulltime, freqs_fd_fulltime, num_freqs_fulltime);
   
    // clear high frequencies
    std::complex<double> y0 =  imag(resp_fd_fulltime[num_freqs_fulltime-1]) + 1i*real(resp_fd_fulltime[num_freqs_fulltime-1]);

    double dy0 =  0;

    double y1 = 0;
    double dy1 = 0;
    double lj;

    double xx = log(spec_len-1)-log(num_freqs_fulltime-1);
    cube_spline0(&sr   ,xx ,real(y0)   ,y1   ,dy0   ,dy1  );
    cube_spline0(&si   ,xx ,imag(y0)   ,0   ,0   ,0  );
   
    std::memset(&fft,0,sizeof(fft));
    init_fft(&fft, spec_len);
    
    // store in Fourier
    for(i = 0; i<num_freqs_fulltime; i++)
        fft.xn[2*i+1] = impulse_spec[2*i+1]*(freqs_fd_fulltime[i]*(imag(resp_fd_fulltime[i]) + 1i*real(resp_fd_fulltime[i])))*2.*M_PI*100.;

   for(i=2*num_freqs_fulltime; i<spec_len; i+=2) 
   {
      lj = log(i)-log(2*num_freqs_fulltime);
      fft.xn[i] = sr.a[3]*lj*lj*lj+sr.a[2]*lj*lj+sr.a[1]*lj+sr.a[0] +
                            1i * (si.a[3]*lj*lj*lj+si.a[2]*lj*lj+si.a[1]*lj+si.a[0]);                     
   }
       
   fft_pro(&fft,1);
   
   // cut into channels
   for(i = 0; i < num_channels; i++)
   {
      dest[i] = 0;
      for(j=chans[i]; j < chans[i+1]; j++)
         dest[i] += real(fft.fn[base_chan + j]);
            
      dest[i] = dest[i]/(chans[i+1] - chans[i]);
   }
   
}

// calculate forward response
void aem_model_raw::fd_forward(std::complex<double> *dest, double *freqs, int num_freqs)
{
    using namespace std::complex_literals;
    int i, j;
    double rho_inf, m, tau, c;
    std::complex<double> alp;
    for(i=0; i<num_freqs; i++)
    {
        for(j=0; j<num_layers; j++) rhos_fd[j] = rhos[j];

        // calculate f dependent resistivities
        for(j=0; j<num_pol_layers; j++)
        {
            rho_inf = cole[3*j];
            m = 1 - rho_inf/rhos[pol_inds[j]];
            tau = cole[3*j + 1];
            c = cole[3*j + 2];
            alp = pow(freqs[i]*tau*2*M_PI, c);
            alp = 1. + alp * std::exp(c*1i*M_PI/2.);
            rhos_fd[pol_inds[j]] = rhos_fd[pol_inds[j]]*(1. - m*(1. - 1./alp));
        }

        dest[i] = ImHz(hor_dist, 2*(alt + altitude_correction) + ver_dist, freqs[i])/prim_field;
    }
}

std::complex<double> aem_model_raw::ImHz(double r,double z,double f)
{
    return integral(z,r,f);
}

std::complex<double> aem_model_raw::integral(double hh, double r, double f)
{
    using namespace std::complex_literals;
    std::complex<double> PS;     
    std::complex<double> intl = 0;
    double dn0;
    double n0=0;
    std::complex<double> n1,c;
    std::complex<double> sigma = 1./rhos_fd[0];
    std::complex<double> Imp;
    double om = f*2*M_PI;
    c = 1i*om*sigma*mu0;   

    double VAL = .001;
    for(n0=VAL,dn0=VAL;n0<1;n0+=dn0) 
    {
        n1 = std::sqrt(n0*n0+c);
        Imp = Impedance(n0,om);
        PS  = PartSum(n0,hh,r,n1,Imp);
        if(std::isnan(real(PS)))
        {
           if(std::isnan(imag(PS)))
             PS = 0. + 0.*1i;
           else
             PS = 0. + imag(PS)*1i;
        } 
        else 
        {
            if(std::isnan(imag(PS)))
              PS = real(PS) + 0.*1i;
 
        }
        intl += dn0*PS;
    }
    return intl;
}

std::complex<double> aem_model_raw::Impedance(double n0, double om)
{
    using namespace std::complex_literals;
    int i, m;
    std::complex<double> ni,nim1;
    std::complex<double> Imp;
    double dpth;

    Imp = 1;
    nim1 = std::sqrt(n0*n0 + 1i*om*mu0/rhos_fd[num_layers-1]);
    dpth = 0;
    m = num_layers - 1;

    for(i=m;i>0;i--)
    {   
       dpth+=depths[i-1];
       ni = nim1;
       nim1 = std::sqrt(n0*n0 + 1i*om*mu0/rhos_fd[i-1]);
       Imp = std::tanh(nim1*dpth+std::atanh(nim1/ni*Imp));
       dpth = 0;
    }

    return Imp;
}

std::complex<double> aem_model_raw::PartSum(double n0, double hh, double r, std::complex<double> n1, std::complex<double> Imp)
{
    std::complex<double> s;
    
    if(fabs(r)>.001)
        s = bessj0(n0*r);
    else 
        s  = 1;

    std::complex<double> A = exp(-n0*hh)*(n1-n0*Imp)*n0*n0*.25/(n1+n0*Imp)/M_PI;
    s = A*s;
    return s;
}

void aem_model_raw::print_model()
{
    int i;

    std::cout << "\n-----------------------------------------------------------------------------------------\n";
    std::cout << "\nRESISTIVITY\n";
    for(i=0; i<num_layers; i++)
        std::cout << rhos[i] << " ";

    std::cout << "\nAC\n";
    std::cout << altitude_correction;

    std::cout << "\nDEPTH\n";
    for(i=0; i<num_layers-1; i++)
        std::cout << depths[i] << " ";

    std::cout << "\nCOLE\n";
    for(i=0; i<num_pol_layers*3; i++)
        std::cout << cole[i] << " ";
    std::cout << "\n-----------------------------------------------------------------------------------------\n";
}

void aem_model_raw::print_to_file(std::ofstream &out)
{
    int i;
    double depp;


    for(i=0; i<num_layers; i++)
        out << rhos[i] << " ";

    out << altitude_correction << " ";

    depp = 0;
    for(i=0; i<num_layers-1; i++)
    {
        out << depp+depths[i]*.5 << " ";
        depp+=depths[i];
    }
}