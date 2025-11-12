#include "spline.h"

void cube_spline(spline *s,double x0, double x1,double y0,double y1,double dy0,double dy1)
{
     s->a[3] = ( 2*(x0*dy0-y0-x1*dy1+y1)+(dy0-dy1)*(x0+x1) )/((x0-x1)*(7*x0*x0+10*x0*x1+7*x1*x1) );
     s->a[2] = ( dy0-dy1-3*(x0-x1)*(x0+x1)*s->a[3] )/( 2.*(x0-x1) );
     s->a[1] = ( (y0-y1) - (x0*x0*x0-x1*x1*x1)*s->a[3] - (x0-x1)*(x0+x1)*s->a[2] )/( x0-x1 );
     s->a[0] = y0 - x0*x0*x0*s->a[3] - x0*x0*s->a[2] - x0*s->a[1];
}

void cube_spline0(spline *s,double x, double y0,double y1,double dy0,double dy1)
{
     s->a[0] = y0;
     s->a[1] = dy0;
     s->a[2] = ((3*(y1-y0)/x)  - (2*dy0+dy1))/x;
     s->a[3] = ((dy0+dy1) - (2*(y1-y0)/x))/x/x;
}