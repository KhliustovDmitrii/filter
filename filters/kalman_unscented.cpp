#include <cmath>
#include "kalman_unscented.h"

// inner structure initialization
kalman_unscented::kalman_unscented(model<double, double> *mod)
{
   m = mod;
   P_y = new double[m->forward_size*m->forward_size];
   P_aux = new double[m->forward_size*m->forward_size];
   P_y_aux = new double[m->forward_size*m->forward_size];
   P_xy = new double[m->num_pars*m->forward_size];
   K = new double[m->num_pars*m->forward_size];
   P = new double[m->num_pars*m->num_pars];
   S_aux_1 = new double[m->num_pars*m->forward_size];
   
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
   // compute sigma points with rows of P
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
      m->conv_to_raw();
      m->response(&(y_sigma[i*m->forward_size]));
   }
   
   // restore the state
   memcpy(m->params, model_state, m->num_pars*sizeof(double));
   m->conv_to_raw();
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
         
   invert_matrix(P_aux, P_y_aux, m->forward_size);
   
   for(i=0; i<m->num_pars; i++)
      for(j=0; j<m->forward_size; j++)
      {
         K[i*m->forward_size + j] = 0;
         for(k=0; k<m->forward_size; k++)
            K[i*m->forward_size + j] = K[i*m->forward_size + j] + P_xy[i*m->forward_size + k]*P_aux[k*m->forward_size + j];
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
   memset(S_aux_1, 0, m->num_pars*m->num_pars*sizeof(double));
   
   // K P_y
   for(i=0; i<m->num_pars; i++)
      for(j=0; j<m->forward_size; j++)
      {
         S_aux_1[i*m->forward_size + j] = 0;
         for(k=0; k<m->forward_size; k++)
            S_aux_1[i*m->forward_size + j] = S_aux_1[i*m->forward_size + j] + K[i*m->forward_size + k]*P_y[k*m->forward_size + j];
      }
     
   // K P_y Kt
   for(i=0; i<m->num_pars; i++)
      for(j=0; j<m->num_pars; j++)
      {
         upd_cov[i*m->num_pars + j] = 0;
         for(k=0; k<m->forward_size; k++)
            upd_cov[i*m->num_pars + j] = upd_cov[i*m->num_pars + j] + S_aux_1[i*m->forward_size + k]*K[j*m->forward_size + k];
      }
}

/* update covariances */
void kalman_unscented::update_cov(double *upd_cov, double factor)
{
   int i, count, j, k;
   
   
   memset(P, 0, m->num_pars*m->num_pars*sizeof(double));
   // compute P = St S
   for(i=0; i<m->num_pars; i++)
      for(j=0; j<m->num_pars; j++)
         for(k=0; k<m->num_pars; k++)
            P[i*m->num_pars + j] = P[i*m->num_pars + j] + S[k*m->num_pars + i]*S[k*m->num_pars + j];
   
   // P - K P_y Kt
   for(i=0; i<m->num_pars; i++)
      for(j=0; j<m->num_pars; j++)
         P[i*m->num_pars + j] = P[i*m->num_pars + j] - upd_cov[i*m->num_pars + j]*factor*factor;
            
   // go back to square root of P = St S
   matrix_square_root(S, P, m->num_pars);
}

void kalman_unscented::invert_matrix(double *dest, double *source, int dim)
{
   int i, j, k;
   double tmp;
   
   for(i=0; i<dim; i++)
      for(j=0; j<dim; j++)
         dest[i*dim + j] = (i==j)?1:0;
         
   for(i=0; i<dim; i++) // use all rows
   {
      // current diagonal element is zero
      // change with row with non-zero
      j=i;
      while(source[j*dim+i]==0) j++;
      
      for(k=0; k<dim; k++)
      {
         tmp = source[j*dim + k];
         source[j*dim + k] = source[i*dim + k];
         source[i*dim + k] = tmp;
         
         tmp = dest[j*dim + k];
         dest[j*dim + k] = dest[i*dim + k];
         dest[i*dim + k] = tmp;
      }
      
      
      // nonzero diag - main part
      for(j=0; j<dim; j++) // all rows but current
      {
         if(i==j) continue;
         
         if(source[j*dim+i] == 0) continue;
         
         tmp = - source[j*dim+i]/source[i*dim+i];
         
         for(k=i; k<dim; k++)
            source[j*dim+k] = source[j*dim+k] + tmp*source[i*dim+k];
            
         for(k=0; k<dim; k++)
            dest[j*dim+k] = dest[j*dim+k] + tmp*dest[i*dim+k];
         
      }
      
      // normalize current row
      tmp = 1/source[i*dim+i];
      for(k=i; k<dim; k++)
         source[i*dim+k] = source[i*dim+k]*tmp;
         
      for(k=0; k<dim; k++)
         dest[i*dim+k] = dest[i*dim+k]*tmp;
   }
   
}

// NOTE: do we need to transpose here?
void kalman_unscented::matrix_square_root(double *dest, double *source, int dim)
{
   int i, j, k;
   
   for(i=0; i<dim; i++)
   {
      dest[i*dim + i] = source[i*dim + i];
      for(j=0; j<i; j++)
         dest[i*dim + i] = dest[i*dim + i] - dest[i*dim + j]*dest[i*dim + j];
      dest[i*dim + i] = std::sqrt(dest[i*dim + i]);
      
      for(j=0; j<i; j++)
         dest[j*dim + i] = 0;
         
      for(j=i+1; j<dim; j++)
      {
         dest[j*dim + i] = source[j*dim + i];
         for(k=0; k<i; k++)
            dest[j*dim + i] = dest[j*dim + i] - dest[j*dim + k]*dest[i*dim + k];
            
         dest[j*dim + i] = dest[j*dim + i]/dest[i*dim + i];
      }
   } 
}
