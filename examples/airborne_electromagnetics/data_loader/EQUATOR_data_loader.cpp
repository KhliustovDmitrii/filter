#include <algorithm>

#include "EQUATOR_data_loader.h"
#include "utilities/mathematics/special_functions/special_functions.h"

void EQUATOR_data_loader::load(const std::vector<double> & data)
{
    int i;

    // load model parameters
    E->E->hor_dist = data[hor_dist_pos];
    E->E->ver_dist = data[ver_dist_pos];
    E->E->alt = data[alt_pos];

    // load measurements

    // fd re
    std::copy(data.begin() + freq_re_pos, 
              data.begin() + freq_re_pos + E->E->freqs.size(),
              measurements.begin());

    // fd im
    std::copy(data.begin() + freq_im_pos, 
              data.begin() + freq_im_pos + E->E->freqs.size(),
              measurements.begin() + E->E->freqs.size());

    // td
    std::copy(data.begin() + chan_pos, 
              data.begin() + chan_pos + E->E->chans.size()-1,
              measurements.begin() + 2*E->E->freqs.size());

    // account for signs
    for(i=0; i<E->E->freqs.size(); i++)
        measurements[i] *= sign_re;
      
    for(i=E->E->freqs.size(); i<2*E->E->freqs.size(); i++)
        measurements[i] *= sign_im;

    for(i=2*E->E->freqs.size(); i<2*E->E->freqs.size() + E->E->chans.size()-1; i++)
        measurements[i] *= sign_td;

    // take difference of in-phase component
    for(i=0; i<E->E->num_freqs-1; i++)
        measurements[i] = measurements[i+1] - measurements[i];
    
    // for highest freq duplicate the previous value
    if(E->E->num_freqs > 1)
        measurements[i] = measurements[i-1];

    // convert to log-linear scale
    for(i=0; i<measurements.size(); i++)
        measurements[i] = log_lin(measurements[i]);
}