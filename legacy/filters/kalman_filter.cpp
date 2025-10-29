#include <cmath>
#include "kalman_filter.h"

kalman_filter::kalman_filter(model<double, double> *mod)
{
   m = mod;
   Jacobian = new double[m->forward_size*m->num_pars];
}

void kalman_filter::get_jacobian()
{
   int i, j;
   double h = 0.001;
   double resp_ini[m->forward_size], resp_var[m->forward_size];
   
   m->response(resp_ini);
   
   for(i=0; i<m->num_pars; i++)
   {
      m->params[i] += h;
      m->conv_to_raw();
      
      m->response(resp_var);
      for(j=0; j<m->forward_size; j++)
         Jacobian[j*m->num_pars + i] = (resp_var[j] - resp_ini[j])/h;
      
      m->params[i] -= h;
      m->conv_to_raw();
   }
}

/*
void kalman_filter::update_model(double *upd_vec, double *resp, double res)
{
   int i, count;
   double resp_tmp[m->forward_size];
   double old_pars[m->num_pars];
   double res_new;
   
   memcpy(old_pars, m->params, sizeof(double)*m->num_pars);
   
   // try to reduce the residual
   for(count=0; count<10; count++)
   {
      // update model params
      for(i=0; i<m->num_pars; i++)
      {
         m->params[i] = old_pars[i] - upd_vec[i];
         
         if(std::isnan(m->params[i])) m->params[i] = 0;
         if(m->params[i] > m->upper_boundary[i]) m->params[i] = m->upper_boundary[i];
         if(m->params[i] < m->lower_boundary[i]) m->params[i] = m->lower_boundary[i];
      }
      
      // calculate new residual
      m->response(resp_tmp);
      res_new = 0;
      for(i=0; i<m->forward_size; i++)
         res_new+=(resp[i]-resp_tmp[i])*(resp[i]-resp_tmp[i])/R[i*m->num_pars + i];
      
      // try to reduce the step
      if(res_new > res*1.1)
         for(i=0; i<m->num_pars; i++)
            upd_vec[i] = upd_vec[i]*0.7;
      else
         break;
   }
   
}
*/
void kalman_filter::get_update_vector(double *upd_vec, double *measurements)
{
   int i;
   double resp[m->forward_size];
   
   kalman_filter::get_jacobian();
   m->response(resp);
   
   std::memset(upd_vec, 0, sizeof(double)*m->num_pars);
   for(i=0; i<m->forward_size; i++)
      kalman_filter::proc(upd_vec, resp[i] - measurements[i], &Jacobian[i*m->num_pars], R[i*m->forward_size + i]*R[i*m->forward_size + i]);
}

void kalman_filter::proc(double *upd_vec, double delta, double *h, double r)
{
   int i,j,k;
   double d[2], bk,ck,dz, tmp;
    
   double *f, *e;
   f = new double[m->num_pars+1];
   e = new double[m->num_pars+1];
    
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
    
   delete(f);
   delete(e);
}




