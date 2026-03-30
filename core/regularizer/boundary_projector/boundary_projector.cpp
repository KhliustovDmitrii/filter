#include <cmath>
#include "boundary_projector.h"

namespace filter::components
{
void Boundary_Projector::prune()
{
    int i;
    double param;

    for(i=0; i<m.num_pars; i++)
    {
        param = m.get_param(i);

        if(std::isnan(param)) param = 0;
        if(param > upper_bound[i]) param = upper_bound[i];
        if(param < lower_bound[i]) param = lower_bound[i];

        m.set_param(i, param);
    }
}
}; // filter::components