#include "kalman_extended.h"

void Kalman_Extended::get_update(double *mes, double *upd_vec, double *upd_cov) 
{
    int i;
    double resp[m->forward_size];
   
    Kalman_Extended::get_jacobian();
    m->response(resp);
   
    std::memset(upd_vec, 0, sizeof(double)*m->num_pars);
    std::memcpy(upd_cov, 0, sizeof(double)*m->num_pars*m->num_pars);

    // process each observation and update parameters accordingly
    for(i=0; i<m->forward_size; i++)
        Kalman_Extended::proc(i, upd_vec, upd_cov, resp[i] - mes[i]);
}

void Kalman_Extended::get_jacobian()
{
    int i, j;
    double h = 0.001;
    double param;

    double resp_ini[m->forward_size], resp_var[m->forward_size];
   
    m->response(resp_ini);
   
    for(i=0; i<m->num_pars; i++)
    {
        // vary the parameter
        param = m->get_param(i);
        m->set_param(i, param+h);

        m->response(resp_var);
        for(j=0; j<m->forward_size; j++)
            Jacobian[j*m->num_pars + i] = (resp_var[j] - resp_ini[j])/h;
      
        // restore original value
        m->set_param(i, param);
    }
}

void Kalman_Extended::proc(int i, double *upd_vec, double *upd_cov, double delta)
{
    int i,j,k;
    double d[2], bk,ck,dz, tmp;

    double f[m->num_pars];
    double e[m->num_pars];

    double *h = &Jacobian[i*m->num_pars];
    double r = R[i*m->forward_size + i]*R[i*m->forward_size + i];
    
    std::memset(e, 0, sizeof(double)*m->num_pars);

    for(i=0; i<m->num_pars; i++) 
    {
        f[i] = 0;
        for(j=0;j<m->num_pars;j++)
            f[i]+=upd_cov[j*m->num_pars+i]*h[j];
    }
   
    d[0] = r;
    for(i=0; i<m->num_pars; i++) 
    {
        d[1] = d[0] + f[i]*f[i];
        bk = sqrt(d[0]/d[1]);
        ck = f[i]/sqrt(d[0]*d[1]);
        for(j=0; j<m->num_pars; j++) 
        {
            tmp = upd_cov[j*m->num_pars+i]*f[i];
            upd_cov[j*m->num_pars+i] = bk*upd_cov[j*m->num_pars+i]-ck*e[j];
            e[j] += tmp;
        }
        d[0] = d[1];
    }
   
    dz = delta;
    for(i=0; i<m->num_pars; i++) dz -= h[i]*upd_vec[i];
   
    dz/=d[0];
    for(i=0; i<m->num_pars; i++) upd_vec[i] += e[i]*dz;
}
