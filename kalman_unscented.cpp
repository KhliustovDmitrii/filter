#include <cmath>
#include "kalman_unscented.h"

// inner structure initialization
kalman_unscented::kalman_unscented(coil_model *mod)
{
   m = mod;
   P_y = new double[m->forward_size*m->forward_size];
   P_y_inv = new double[m->forward_size*m->forward_size];
   P_y_aux = new double[m->forward_size*m->forward_size];
   P_xy = new double[m->num_pars*m->forward_size];
   K = new double[m->num_pars*m->forward_size];
   P = new double[m->num_pars*m->num_pars];
   S_aux = new double[m->num_pars*m->forward_size];
   
   x_sigma = new double[2*m->num_pars*m->num_pars];
   y_sigma = new double[2*m->num_pars*m->forward_size];
   
   y_avg = new double[m->forward_size];
   
   P = new double[m->num_pars*m->num_pars];
}

// sigma point computation
void kalman_unscented::compute_sigma()
{
   int i, j;
   double sq_n = std::sqrt(m->num_pars);
   double model_state[m->num_pars];
   
   // compute sigma points with rows of S
   for(i=0; i<m->num_pars; i++)
   {
      for(j=0; j<m->num_pars; j++)
      {
         x_sigma[i*m->num_pars + j] = m->params[j] + sq_n*S[i*m->num_pars + j];
         x_sigma[(i+m->num_pars)*m->num_pars + j] = m->params[j] - sq_n*S[i*m->num_pars + j];
      }
   }
   
   // compute sigma responses
   
   // save model state
   memcpy(model_state, m->params, m->num_pars*sizeof(double));
   
   for(i=0; i<2*m->num_pars; i++)
   {
      memcpy(m->params, &(x_sigma[i*m->num_pars]), m->num_pars*sizeof(double));
      m->response(&(y_sigma[i*m->forward_size]));
   }
   
   // restore the state
   memcpy(m->params, model_state, m->num_pars*sizeof(double));
}

// necessary matrix computations
void kalman_unscented::get_matrices()
{
   int i, j, k;
   
   // compute average response
   memset(y_avg, 0, m->forward_size*sizeof(double));
   for(i=0; i<2*m->num_pars; i++)
      for(j=0; j<m->forward_size; j++)
         y_avg[j] = y_avg[j] + y_sigma[i*m->forward_size + j]/(2*m->num_pars);
   
   // response covariance
   for(i=0; i<m->forward_size; i++)
      for(j=0; j<m->forward_size; j++)
      {
         P_y[i*m->forward_size + j] = R[i*m->forward_size + j];
         for(k=0; k<2*m->num_pars; k++)
            P_y[i*m->forward_size + j] = P_y[i*m->forward_size + j] + 
            (y_sigma[k*m->forward_size + i] - y_avg[i])*(y_sigma[k*m->forward_size + j] - y_avg[j])/(2*m->num_pars);
      }
   
   // cross-covariance
   for(i=0; i<m->num_pars; i++)
      for(j=0; j<m->forward_size; j++)
      {
         P_xy[i*m->forward_size + j] = 0;
         for(k=0; k<2*m->num_pars; k++)
            P_xy[i*m->forward_size + j] = P_xy[i*m->forward_size + j] + 
            (x_sigma[k*m->num_pars + i] - m->params[i])*(y_sigma[k*m->forward_size + j] - y_avg[j])/(2*m->num_pars);
      }
   
   // Kalman gain
   for(i=0; i<m->forward_size; i++)
      for(j=0; j<m->forward_size; j++)
         P_y_aux[i*m->forward_size + j] = P_y[i*m->forward_size + j];
         
   invert_matrix(P_y_inv, P_y_aux, m->forward_size);
   
   for(i=0; i<m->num_pars; i++)
      for(j=0; j<m->forward_size; j++)
      {
         K[i*m->forward_size + j] = 0;
         for(k=0; k<m->forward_size; k++)
            K[i*m->forward_size + j] = K[i*m->forward_size + j] + P_xy[i*m->forward_size + k]*P_y_inv[k*m->forward_size + j];
      }

}

// get updates for model parameters and covariances
void kalman_unscented::get_update_vector(double *upd_vec, double *upd_cov, double *measurements)
{
   int i, j, k;
   
   memset(upd_vec, 0, m->num_pars*sizeof(double));
   memset(upd_cov, 0, m->num_pars*m->num_pars*sizeof(double));
   
   // parameter update
   for(i=0; i<m->num_pars; i++)
      for(j=0; j<m->forward_size; j++)
         upd_vec[i] = upd_vec[i] + K[i*m->forward_size + j]*(measurements[j] - y_avg[j]);
         
   // covariance update
   memset(S_aux, 0, m->num_pars*m->num_pars*sizeof(double));
   
   // K P_y
   for(i=0; i<m->num_pars; i++)
      for(j=0; j<m->forward_size; j++)
      {
         S_aux[i*m->forward_size + j] = 0;
         for(k=0; k<m->forward_size; k++)
            S_aux[i*m->forward_size + j] = S_aux[i*m->forward_size + j] + K[i*m->forward_size + k]*P_y[k*m->forward_size + j];
      }
     
   // K P_y Kt
   for(i=0; i<m->num_pars; i++)
      for(j=0; j<m->num_pars; j++)
      {
         upd_cov[i*m->num_pars + j] = 0;
         for(k=0; k<m->forward_size; k++)
            upd_cov[i*m->num_pars + j] = upd_cov[i*m->num_pars + j] + S_aux[i*m->forward_size + k]*K[j*m->forward_size + k];
      }
}
