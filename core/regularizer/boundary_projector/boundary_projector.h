#ifndef _BOUNDARY_PROJECTOR_H_
#define _BOUNDARY_PROJECTOR_H_

#include <vector>

#include "types/regularizer/regularizer.h"
#include "types/model/model.h"

// detects model parmeters exceeding boundaries and crops them
class Boundary_Projector : public Regularizer
{
private:
    std::vector<double> upper_bound, lower_bound;

public:
    Boundary_Projector(Model *model, std::vector<double> ub, std::vector<double> lb) : 
    Regularizer(model), 
    upper_bound(std::move(ub)), lower_bound(std::move(lb)) {};

    void prune();
};

#endif