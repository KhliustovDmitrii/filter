#ifndef _MATH_UTILS_H_
#define _MATH_UTILS_H_

typedef struct {         // spline structure for interpolation
     double a[4];
} spline;

void cube_spline(spline *s,double x0, double x1,double y0,double y1,double dy0,double dy1);
void cube_spline0(spline *s,double x, double y0,double y1,double dy0,double dy1);
double bessj0(double x);
double log_lin(double x);
double exp_lin(double x);
void var_mes(double *S0, double *S1, double *mes_buf, int sdim);
double primField(double hd, double vd);
#endif