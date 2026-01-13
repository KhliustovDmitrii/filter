#include "decay_updater.h"

double Decay_Updater::update(std::vector<double> &upd_vec, std::vector<double> &mes)
{
    int i, count;
    std::vector<double> resp(m->forward_size, 0);
    std::vector<double> old_pars = std::vector<double>(m->num_pars);
    double res_old, res_new;
    double total_factor = 1;
   
    for(i=0; i<m->num_pars; i++) old_pars[i] = m->get_param(i);

    m->response(resp);
    res_old = m->residual(mes, resp);
   
    // try to reduce the residual
    for(count=0; count<num_steps; count++)
    {
        // update model params
        for(i=0; i<m->num_pars; i++)
            m->set_param(i, old_pars[i] - upd_vec[i]);

        // calculate new residual
        m->response(resp);
        res_new = m->residual(mes, resp);
      
        // try to reduce the step
        if(res_new > res_old)
        {
            for(i=0; i<m->num_pars; i++)
                upd_vec[i] = upd_vec[i]*factor;

            total_factor *= factor;
        }
        else break;
   }

   return total_factor;
}