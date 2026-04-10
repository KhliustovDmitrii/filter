#include "kalman_extended.h"

namespace filter
{

std::unique_ptr<Filter_Workspace> Kalman_Extended::allocate_workspace() const
{
    auto ws = std::make_unique<EKF_Workspace>();

    ws->Jacobian.resize(m.num_pars * m.forward_size, 0.0);
    ws->f.resize(m.num_pars, 0.0);
    ws->e.resize(m.num_pars, 0.0);
    return ws;
}

void Kalman_Extended::get_update(std::vector<double> &mes, 
                                 std::vector<double> &upd_vec, 
                                 std::vector<double> &upd_cov,
                                 Filter_Workspace &ws) const
{
    int i;
    std::vector<double> resp(m.forward_size, 0);

    auto& ekf_ws = static_cast<EKF_Workspace&>(ws);
   
    get_jacobian(ekf_ws);
    m.response(resp);
   
    std::fill(upd_vec.begin(), upd_vec.end(), 0);
    std::copy(S.begin(), S.end(), upd_cov.begin());

    // process each observation and update parameters accordingly
    for(i=0; i<m.forward_size; i++)
        proc(ekf_ws, i, upd_vec, upd_cov, resp[i] - mes[i]);
}

void Kalman_Extended::get_jacobian(EKF_Workspace &ws) const
{
    int i, j;
    double h = 0.001;
    double param;

    std::vector<double> resp_ini(m.forward_size, 0);
    std::vector<double> resp_var(m.forward_size, 0);
   
    m.response(resp_ini);
   
    for(i=0; i<m.num_pars; i++)
    {
        // vary the parameter
        param = m.get_param(i);
        m.set_param(i, param+h);

        m.response(resp_var);
        for(j=0; j<m.forward_size; j++)
            ws.Jacobian[j*m.num_pars + i] = (resp_var[j] - resp_ini[j])/h;
      
        // restore original value
        m.set_param(i, param);
    }
}

void Kalman_Extended::proc(EKF_Workspace &ws, int k, 
                           std::vector<double> &upd_vec, std::vector<double> &upd_cov, double delta) const
{
    int i,j;
    double d[2], bk,ck,dz, tmp;

    double *h = &ws.Jacobian[k*m.num_pars];
    double r = R[k*m.forward_size + k]*R[k*m.forward_size + k];
    std::fill(ws.e.begin(), ws.e.end(), 0);
    std::fill(ws.f.begin(), ws.f.end(), 0);

    upd_cov = S;

    for(i=0; i<m.num_pars; i++) 
    {
        for(j=0;j<m.num_pars;j++)
            ws.f[i]+=upd_cov[j*m.num_pars+i]*h[j];
    }
   
    d[0] = r;
    for(i=0; i<m.num_pars; i++) 
    {
        d[1] = d[0] + ws.f[i]*ws.f[i];
        bk = sqrt(d[0]/d[1]);
        ck = ws.f[i]/sqrt(d[0]*d[1]);
        for(j=0; j<m.num_pars; j++) 
        {
            tmp = upd_cov[j*m.num_pars+i]*ws.f[i];
            upd_cov[j*m.num_pars+i] = bk*upd_cov[j*m.num_pars+i]-ck*ws.e[j];
            ws.e[j] += tmp;
        }
        d[0] = d[1];
    }
   
    dz = delta;
    for(i=0; i<m.num_pars; i++) dz -= h[i]*upd_vec[i];
   
    dz/=d[0];
    for(i=0; i<m.num_pars; i++) upd_vec[i] += ws.e[i]*dz;
}
}; // filter