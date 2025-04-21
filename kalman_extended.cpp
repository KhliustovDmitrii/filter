#include <cmath>
#include "kalman_extended.h"

kalman_extended::kalman_extended(model *mod)
{
   m = mod;
   Jacobian = new double[m->forward_size*m->num_pars];
   f = new double[m->num_pars+1];
   e = new double[m->num_pars+1];
}

void kalman_extended::get_update_vector(double *upd_vec, double *measurements
                                        int *valid_mes)
{
   int i;
   double resp[m->forward_size];
   
   // TODO
   get_jacobian(Jacobian);
   
   m->response(resp);
   
   std::memset(upd_vec, 0, sizeof(double)*m->num_pars);
   for(i=0; i<m->forward_size; i++)
      if(valid_mes[i])
         proc(upd_vec, resp[i] - measurements[i], &Jacobian[i*m->num_pars], R[i*m->forward_size + i]);
}

void kalman_extended::proc(double *upd_vec, double delta, double *h, double r)
{
   int i,j,k;
   double d[2], bk,ck,dz, tmp;
    
   std::memset(e, 0, sizeof(double)*(m->num_pars+1));
   std::memset(f, 0, sizeof(double)*(m->num_pars+1));
    
   d[0] = r;
   for(i=0; i<m->num_pars; i++) 
   {
      f[i] = 0;
      for(j=0;j<m->num_pars;j++)
         f[i]+=S[j*m->num_pars+i]*h[j];
   }
   
   for(i=0; i<m->num_pars; i++) 
   {
      d[1] = d[0] + f[i]*f[i];
      bk = sqrt(d[0]/d[1]);
      ck = f[i]/sqrt(d[0]*d[1]);
      for(j=0; j<m->num_pars; j++) 
      {
         tmp = S[j*m->num_pars+i]*f[i];
         S[j*m->num_pars+i] = bk*S[j*m->num_pars+i]-ck*e[j];
         e[j] += tmp;
      }
      d[0] = d[1];
   }
   
   dz = delta;
   for(i=0; i<m->num_pars; i++) dz -= h[i]*upd_vec[i];
   
   dz/=d[0];
   for(i=0; i<m->num_pars; i++) upd_vec[i] += e[i]*dz;
}




