#include <algorithm>

#include "EQUATOR_data_loader.h"

void EQUATOR_data_loader::load(const std::vector<double> & data)
{
    // load model parameters
    E->E->hor_dist = data[hor_dist_pos];
    E->E->ver_dist = data[ver_dist_pos];
    E->E->altitude_correction[0] = data[alt_pos];

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
}