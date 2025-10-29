#include "kalman_unscented.h"
#include "../../utilities/mathematics/matrix.h"

// sigma point computation
void kalman_unscented::compute_sigma()
{
    int i, j;
    double sq_n = std::sqrt(m->num_pars);
    double model_state[m->num_pars];

    // compute sigma points with rows of covariance matrix
    for(i=0; i<m->num_pars; i++)
    {
        for(j=0; j<m->num_pars; j++)
        {
            x_sigma[i*m->num_pars + j] = m->get_param[j] + sq_n*S[i*m->num_pars + j];
            x_sigma[(i+m->num_pars)*m->num_pars + j] = m->get_param(j) - sq_n*S[i*m->num_pars + j];
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
            m->set_param(j, x_sigma[i*m->num_pars+j]);

        m->response(&(y_sigma[i*m->forward_size]));
    }
   
    // restore the state
    for(i=0; i<m->num_pars; i++)
        m->set_param(i, model_state[i]);
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

void get_update(double *mes, double *upd_vec, double *upd_cov)
{
    int i, j, k;
   
    memset(upd_vec, 0, m->num_pars*sizeof(double));
    memset(upd_cov, 0, m->num_pars*m->num_pars*sizeof(double));
   
    // parameter update
    for(i=0; i<m->num_pars; i++)
        for(j=0; j<m->forward_size; j++)
            upd_vec[i] = upd_vec[i] + K[i*m->forward_size + j]*(mes[j] - y_avg[j]);
         
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