#include "boundary_projector.h"

void boundary_projector::prune()
{
    int i;
    double param;

    for(i=0; i<m->num_pars; i++)
    {
        param = m->get_param(i);

        if(std::isnan(param)) param = 0;
        if(param > upper_boundary[i]) param = upper_boundary[i];
        if(param < lower_boundary[i]) param = lower_boundary[i];

        m->set_param(i, param);
    }
}