#include <cmath>

#include "special_functions.h"

double bessj0(double x)
{
   double ax,z;
   double xx,y,ans,ans1,ans2;

   if ((ax=std::abs(x)) < 8.0) {
      y=x*x;
      ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
         +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
      ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
         +y*(59272.64853+y*(267.8532712+y*1.0))));
      ans=ans1/ans2;
   } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-0.785398164;
      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
         +y*(-0.6911147651e-5+y*(0.7621095161e-6
         -y*0.934935152e-7)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
   }
   return ans;
}

double log_lin(double x)
{
   if(x > 1)
      return 1 + log(x);
   if(x > -1)
      return x;
       
   return -1 - log(-x);
}

double exp_lin(double x)
{
   if(x > 1)
      return exp(x-1);
   if(x > -1)
      return x;
       
   return -exp(-x - 1);
}

double primField(double hd,
                 double vd)
{

    double E[10], R[4], RR[10], Hp[4], k, MR, Ampl, M;
    R[1]=hd; R[2]=0; R[3]=vd; //towed cable = 39.11
    //R[1]=33; R[2]=0; R[3]=21; //towed cable = 39.11
    //M[1]=0; M[2]=0; M[3]=1.e9;
    M = 1;
    E[1]=1; E[2]=0; E[3]=0;
    E[4]=0; E[5]=1; E[6]=0;
    E[7]=0; E[8]=0; E[9]=1;
    RR[1]=R[1]*R[1]; RR[2]=R[1]*R[2]; RR[3]=R[1]*R[3];
    RR[4]=R[2]*R[1]; RR[5]=R[2]*R[2]; RR[6]=R[2]*R[3];
    RR[7]=R[3]*R[1]; RR[8]=R[3]*R[2]; RR[9]=R[3]*R[3];
    MR=R[1]*R[1]+R[2]*R[2]+R[3]*R[3];
    k = M / MR / (sqrt(MR)) / 4 / M_PI;
    Hp[1] = (3*RR[3] / MR-E[3]) * k;
    Hp[2] = (3*RR[6] / MR-E[6]) * k;
    Hp[3] = (3*RR[9] / MR-E[9]) * k;

    Ampl = sqrt(Hp[1] * Hp[1] + Hp[2] * Hp[2] + Hp[3] * Hp[3]);
    return Ampl;
}

double transform_function_dispatch(double x, int t)
{
   switch(t)
   {
      case 0:
         return x;
      case 1:
         return log(x);
      case 2:
         return log_lin(x); 
   }
   
   return x;
}

double inverse_transform_function_dispatch(double x, int t)
{
   switch(t)
   {
      case 0:
         return x;
      case 1:
         return exp(x);
      case 2:
         return exp_lin(x); 
   }

   return x;
}