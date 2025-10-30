#ifndef _EQUATOR_C_H_
#define _EQUATOR_C_H_

// compressed representation of EQUATOR
// complies to model interface
#include "../EQUATOR/EQUATOR.h"
#include "../../../types/model/model.h"

class EQUATOR_C : public model
{
public:

    EQUATOR *E;

    int num_diff_depth;
    int num_diff_rho;
    int num_diff_cole;

    // differentiability - indexes of variable parameters
    int diff_depth[num_diff_depth];
    int diff_rho[num_diff_rho];
    int diff_cole[num_diff_cole];
    int diff_alt_cor;

    // scales; 0 for lin, 1 for log, 2 for log-lin
    int scale_depth[num_diff_depth];
    int scale_rho[num_diff_rho];
    int scale_cole[num_diff_cole];
    int scale_alt_cor;
    int scale_total[num_diff_depth + num_diff_rho + num_diff_cole + diff_alt_cor];

private:

    // weights to compute residual with
    double weights[forward_size];

    // transforms underlying EQUATOR to _C representation
    void compress();
}

#endif