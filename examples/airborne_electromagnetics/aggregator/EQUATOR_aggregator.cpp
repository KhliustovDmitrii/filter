#include <vector>

#include "EQUATOR_aggregator.h"

// function for aggregating resistivities in several models
// Side effect: modifies out
void EQUATOR_aggregator::aggregate(std::vector<double> weights)
{
    int lnum, mnum, lcount;
    double rho, depth_m, depth;

    // process all layers of out model
    depth = 0;
    for(lnum=0; lnum<out_model->E->num_layers; lnum++)
    {
        depth = depth + out_model->E->depths[lnum];
        rho = 0;
        for(mnum=0; mnum<source_models.size(); mnum++)
        {
            lcount = 0;
            depth_m = 0;

            // find the layer in source models which covers the current layer
            // in out model
            while(lcount < source_models[mnum]->E->num_layers - 1 &&
                depth_m < depth)
            {
                lcount++;
                depth_m += source_models[mnum]->E->depths[lcount];
            }

            rho = rho + weights[mnum]*source_models[mnum]->E->rhos[lcount];
        }

        out_model->E->rhos[lnum] = rho;
    }

    out_model->full_to_c();
}