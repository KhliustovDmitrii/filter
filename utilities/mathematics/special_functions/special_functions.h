#ifndef SPECIAL_FUNCTIONS_H
#define SPECIAL_FUNCTIONS_H

namespace filter::math
{
double bessj0(double x);
double log_lin(double x);
double exp_lin(double x);
double primField(double hd, double vd);
double transform_function_dispatch(double x, int t);
double inverse_transform_function_dispatch(double x, int t);
}; // filter::math
#endif