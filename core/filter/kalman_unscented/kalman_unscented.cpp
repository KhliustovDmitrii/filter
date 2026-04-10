#include <cmath>

#include "kalman_unscented.h"
#include "utilities/mathematics/matrix/matrix.h"

namespace filter
{

std::unique_ptr<Filter_Workspace> Kalman_Unscented::allocate_workspace() const
{
    auto ws = std::make_unique<UKF_Workspace>();

    ws->P_y.resize(m.forward_size*m.forward_size, 0);
    ws->P_xy.resize(m.num_pars*m.forward_size, 0);
    ws->P_aux.resize(m.forward_size*m.forward_size, 0);
    ws->P.resize(m.num_pars*m.num_pars, 0);
    ws->S_aux_1.resize(m.num_pars*m.forward_size, 0);
    ws->P_y_aux.resize(m.forward_size*m.forward_size);
    ws->K.resize(m.num_pars*m.forward_size, 0);
    ws->x_sigma.resize(2*m.num_pars, std::vector<double>(m.num_pars, 0));
    ws->y_sigma.resize(2*m.num_pars, std::vector<double>(m.forward_size, 0));
    ws->y_avg.resize(m.forward_size, 0);

    return ws;
}

// sigma point computation
void Kalman_Unscented::compute_sigma(Filter_Workspace &ws) const
{
    int i, j;
    double sq_n = sqrt(m.num_pars);
    std::vector<double> model_state(m.num_pars, 0);

    auto& ukf_ws = static_cast<UKF_Workspace&>(ws);

    // compute sigma points with rows of covariance matrix
    for(i=0; i<m.num_pars; i++)
    {
        for(j=0; j<m.num_pars; j++)
        {
            ukf_ws.x_sigma[i][j] = m.get_param(j) + sq_n*S[i*m.num_pars + j];
            ukf_ws.x_sigma[i+m.num_pars][j] = m.get_param(j) - sq_n*S[i*m.num_pars + j];
        }
    }
   
    // compute sigma responses
   
    // save model state
    for(i=0; i<m.num_pars; i++)
        model_state[i] = m.get_param(i);
   
    // compute response
    for(i=0; i<2*m.num_pars; i++)
    {
        for(j=0; j<m.num_pars; j++)
            m.set_param(j, ukf_ws.x_sigma[i][j]);

        m.response(ukf_ws.y_sigma[i]);
    }
   
    // restore the state
    for(i=0; i<m.num_pars; i++)
        m.set_param(i, model_state[i]);
}

// necessary matrix computations
void Kalman_Unscented::get_matrices(Filter_Workspace &ws) const
{
    int i, j, k;
    auto& ukf_ws = static_cast<UKF_Workspace&>(ws);
   
    // compute average response
    std::fill(ukf_ws.y_avg.begin(), ukf_ws.y_avg.end(), 0);
    for(i=0; i<2*m.num_pars; i++)
        for(j=0; j<m.forward_size; j++)
            ukf_ws.y_avg[j] = ukf_ws.y_avg[j] + ukf_ws.y_sigma[i][j]/(2*m.num_pars);
   
    // response covariance
    for(i=0; i<m.forward_size; i++)
        for(j=0; j<m.forward_size; j++)
        {
            ukf_ws.P_y[i*m.forward_size + j] = R[i*m.forward_size + j];
            for(k=0; k<2*m.num_pars; k++)
                ukf_ws.P_y[i*m.forward_size + j] = ukf_ws.P_y[i*m.forward_size + j] + 
                (ukf_ws.y_sigma[k][i] - ukf_ws.y_avg[i])*(ukf_ws.y_sigma[k][j] - ukf_ws.y_avg[j])/(2*m.num_pars);
        }
   
    // cross-covariance
    for(i=0; i<m.num_pars; i++)
        for(j=0; j<m.forward_size; j++)
        {
            ukf_ws.P_xy[i*m.forward_size + j] = 0;
            for(k=0; k<2*m.num_pars; k++)
                ukf_ws.P_xy[i*m.forward_size + j] = ukf_ws.P_xy[i*m.forward_size + j] + 
                (ukf_ws.x_sigma[k][i] - m.get_param(i))*(ukf_ws.y_sigma[k][j] - ukf_ws.y_avg[j])/(2*m.num_pars);
        }
   
    // Kalman gain
    for(i=0; i<m.forward_size; i++)
        for(j=0; j<m.forward_size; j++)
            ukf_ws.P_y_aux[i*m.forward_size + j] = ukf_ws.P_y[i*m.forward_size + j];
         
    math::invert_matrix(ukf_ws.P_aux, ukf_ws.P_y_aux, m.forward_size);
   
    for(i=0; i<m.num_pars; i++)
        for(j=0; j<m.forward_size; j++)
        {
            ukf_ws.K[i*m.forward_size + j] = 0;
            for(k=0; k<m.forward_size; k++)
             ukf_ws.K[i*m.forward_size + j] = ukf_ws.K[i*m.forward_size + j] +ukf_ws.P_xy[i*m.forward_size + k]*ukf_ws.P_aux[k*m.forward_size + j];
        }
}

void Kalman_Unscented::get_update(std::vector<double> &mes, 
    std::vector<double> &upd_vec, 
    std::vector<double> &upd_cov,
    Filter_Workspace &ws) const
{
    int i, j, k;
    auto& ukf_ws = static_cast<UKF_Workspace&>(ws);
   
    std::fill(upd_vec.begin(), upd_vec.end(), 0);
    std::fill(upd_cov.begin(), upd_cov.end(), 0);

    compute_sigma(ws);
    get_matrices(ws);
   
    // parameter update
    for(i=0; i<m.num_pars; i++)
        for(j=0; j<m.forward_size; j++)
            upd_vec[i] = upd_vec[i] - ukf_ws.K[i*m.forward_size + j]*(mes[j] - ukf_ws.y_avg[j]);
         
    // covariance update
    std::fill(ukf_ws.S_aux_1.begin(), ukf_ws.S_aux_1.end(), 0);

    // K P_y
    for(i=0; i<m.num_pars; i++)
        for(j=0; j<m.forward_size; j++)
        {
            ukf_ws.S_aux_1[i*m.forward_size + j] = 0;
            for(k=0; k<m.forward_size; k++)
                ukf_ws.S_aux_1[i*m.forward_size + j] = ukf_ws.S_aux_1[i*m.forward_size + j] + ukf_ws.K[i*m.forward_size + k]*ukf_ws.P_y[k*m.forward_size + j];
        }
     
    // K P_y Kt
    for(i=0; i<m.num_pars; i++)
        for(j=0; j<m.num_pars; j++)
        {
            upd_cov[i*m.num_pars + j] = 0;
            for(k=0; k<m.forward_size; k++)
                upd_cov[i*m.num_pars + j] = upd_cov[i*m.num_pars + j] + ukf_ws.S_aux_1[i*m.forward_size + k]*ukf_ws.K[j*m.forward_size + k];
        }
}

void Kalman_Unscented::update_covariance(Filter_Workspace &ws, std::vector<double> &upd_cov, double factor)
{
   int i, count, j, k;
   auto& ukf_ws = static_cast<UKF_Workspace&>(ws);
   
   std::fill(ukf_ws.P.begin(), ukf_ws.P.end(), 0);

   // compute P = St S
   for(i=0; i<m.num_pars; i++)
      for(j=0; j<m.num_pars; j++)
         for(k=0; k<m.num_pars; k++)
            ukf_ws.P[i*m.num_pars + j] = ukf_ws.P[i*m.num_pars + j] + S[k*m.num_pars + i]*S[k*m.num_pars + j];
   
   // P - K P_y Kt
   for(i=0; i<m.num_pars; i++)
      for(j=0; j<m.num_pars; j++)
         ukf_ws.P[i*m.num_pars + j] = ukf_ws.P[i*m.num_pars + j] - upd_cov[i*m.num_pars + j]*factor*factor;
            
   // go back to square root of P = St S
   math::matrix_square_root(S, ukf_ws.P, m.num_pars);
}
}; // filter