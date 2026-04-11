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

struct LogLin
{
    double operator()(double x)
    {
        if(x > 1)
            return 1 + log(x);
        if(x > -1)
            return x;
            
        return -1 - log(-x);
    }
};

struct ExpLin
{
    double operator()(double x)
    {
        if(x > 1)
            return exp(x-1);
        if(x > -1)
            return x;
            
        return -exp(-x - 1);
    }
};

struct Lin
{
    double operator()(double x)
    {
        return x;
    }
};

}; // filter::math
#endif