#ifndef _SPLINE_H_
#define _SPLINE_H_
typedef struct {         // spline structure for interpolation
     double a[4];
} spline;
void cube_spline(spline *s,double x0, double x1,double y0,double y1,double dy0,double dy1);
void cube_spline0(spline *s,double x, double y0,double y1,double dy0,double dy1);
#endif