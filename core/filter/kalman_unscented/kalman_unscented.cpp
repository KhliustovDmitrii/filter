#include <cmath>

#include "kalman_unscented.h"
#include "utilities/mathematics/matrix/matrix.h"

// sigma point computation
void Kalman_Unscented::compute_sigma()
{
    int i, j;
    double sq_n = sqrt(m->num_pars);
    double model_state[m->num_pars];

    // compute sigma points with rows of covariance matrix
    for(i=0; i<m->num_pars; i++)
    {
        for(j=0; j<m->num_pars; j++)
        {
            x_sigma[i][j] = m->get_param(j) + sq_n*S[i*m->num_pars + j];
            x_sigma[i+m->num_pars][j] = m->get_param(j) - sq_n*S[i*m->num_pars + j];
        }
    }
   
    // compute sigma responses
   
    // save model state
    for(i=0; i<m->num_pars; i++)
        model_state[i] = m->get_param(i);
   
    // compute response
    for(i=0; i<2*m->num_pars; i++)
    {
        for(j=0; j<m->num_pars; j++)
            m->set_param(j, x_sigma[i][j]);

        m->response(y_sigma[i]);
    }
   
    // restore the state
    for(i=0; i<m->num_pars; i++)
        m->set_param(i, model_state[i]);
}

// necessary matrix computations
void Kalman_Unscented::get_matrices()
{
    int i, j, k;
   
    // compute average response
    std::fill(y_avg.begin(), y_avg.end(), 0);
    for(i=0; i<2*m->num_pars; i++)
        for(j=0; j<m->forward_size; j++)
            y_avg[j] = y_avg[j] + y_sigma[i][j]/(2*m->num_pars);
   
    // response covariance
    for(i=0; i<m->forward_size; i++)
        for(j=0; j<m->forward_size; j++)
        {
            P_y[i*m->forward_size + j] = R[i*m->forward_size + j];
            for(k=0; k<2*m->num_pars; k++)
                P_y[i*m->forward_size + j] = P_y[i*m->forward_size + j] + 
                (y_sigma[k][i] - y_avg[i])*(y_sigma[k][j] - y_avg[j])/(2*m->num_pars);
        }
   
    // cross-covariance
    for(i=0; i<m->num_pars; i++)
        for(j=0; j<m->forward_size; j++)
        {
            P_xy[i*m->forward_size + j] = 0;
            for(k=0; k<2*m->num_pars; k++)
                P_xy[i*m->forward_size + j] = P_xy[i*m->forward_size + j] + 
                (x_sigma[k][i] - m->get_param(i))*(y_sigma[k][j] - y_avg[j])/(2*m->num_pars);
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

void Kalman_Unscented::get_update(std::vector<double> &mes, std::vector<double> &upd_vec, std::vector<double> &upd_cov)
{
    int i, j, k;
   
    std::fill(upd_vec.begin(), upd_vec.end(), 0);
    std::fill(upd_cov.begin(), upd_cov.end(), 0);

    compute_sigma();
    get_matrices();
   
    // parameter update
    for(i=0; i<m->num_pars; i++)
        for(j=0; j<m->forward_size; j++)
            upd_vec[i] = upd_vec[i] - K[i*m->forward_size + j]*(mes[j] - y_avg[j]);
         
    // covariance update
    std::fill(S_aux_1.begin(), S_aux_1.end(), 0);

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

void Kalman_Unscented::update_covariance(std::vector<double> &upd_cov, double factor)
{
   int i, count, j, k;
   
   std::fill(P.begin(), P.end(), 0);

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