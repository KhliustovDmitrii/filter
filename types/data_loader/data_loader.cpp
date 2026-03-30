#include "data_loader.h"

namespace filter::io
{
void Data_Loader::get_measurements(std::vector<double> &mes)
{
    std::copy(measurements.begin(), measurements.end(), mes.begin());
}
}; // filter::io