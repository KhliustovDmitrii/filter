#include "data_loader.h"

void Data_Loader::get_measurements(std::vector<double> &mes)
{
    std::copy(measurements.begin(), measurements.end(), mes.begin());
}