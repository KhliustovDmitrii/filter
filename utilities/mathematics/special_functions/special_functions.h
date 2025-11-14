#ifndef _SPECIAL_FUNCTIONS_H_
#define _SPECIAL_FUNCTIONS_H_
double bessj0(double x);
double log_lin(double x);
double exp_lin(double x);
double primField(double hd, double vd);
double transform_function_dispatch(double x, int t);
double inverse_transform_function_dispatch(double x, int t);
#endif