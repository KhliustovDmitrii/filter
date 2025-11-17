#include <algorithm>

#include "./averaging_filter.h"

int Averaging_Filter::add_data(std::vector<double> &data, std::vector<double> &output)
{
    int i;

    for(i=0; i<data.size(); i++)
        mes_state[i] += data[i];

    agg_state++;

    // time to average
    if(agg_state == window_size)
    {
        for(i=0; i<data.size(); i++)
            output[i] = mes_state[i]/window_size;

        agg_state = 0;
        std::fill(mes_state.begin(), mes_state.end(), 0);
        return 0;
    }

    return 1;
}