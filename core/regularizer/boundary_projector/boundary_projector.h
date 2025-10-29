#ifndef _BOUNDARY_PROJECTOR_H_
#define _BOUNDARY_PROJECTOR_H_

#include "../types/regularizer/regularizer.h"
#include "../types/model/model.h"

class boundary_projector : public regularizer
{
private:
    double upper_bound[m->num_pars];
    double lower_bound[m->num_pars];
}

#endif