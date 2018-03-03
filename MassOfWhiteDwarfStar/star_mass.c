#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cpgplot.h"


double f(double x, double y, double t)
{
 /* DE for density */
 double vgamma(double x)
 {
  x = pow(fabs(x),1.0/3.0);
  return x*x/(3*sqrt(1+x*x));
 }
 return -x*y/(vgamma(x)*t*t);
}

double g(double x, double y, double t)
{
 /* DE for mass */
 return t*t*x;
}

int rk_coupled(double * coord, double step, double (* f) (double,double,double), double (* g) (double,double,double))
{
 /* Updates the coordinates by 1 step. Important note: FOr the first parameter,
    need to point coord to the first coordinate first */
 double xn = *coord;
 coord++;
 double yn = *coord;
 coord++;
 double tn = *coord;
 /* Finding the gradients */
 double f1 = f(xn,yn,tn);
 double g1 = g(xn,yn,tn);
 double f2 = f((xn + (step/2)*f1), yn + ((step/2)*g1), tn + step/2);
 double g2 = g((xn + (step/2)*f1), yn + ((step/2)*g1), tn + step/2);
 double f3 = f((xn + (step/2)*f2), yn + ((step/2)*g2), tn + step/2);
 double g3 = g((xn + (step/2)*f2), yn + ((step/2)*g2), tn + step/2);
 double f4 = f(xn + step*f3, yn + step*g3, tn + step);
 double g4 = g(xn + step*f3, yn + step*g3, tn + step);
 * coord = tn + step;
 coord = coord - 1;
 * coord = yn + (step/6)*(g1 + 2*g2 + 2*g3 + g4);
 coord = coord - 1;
 * coord = xn + (step/6)*(f1 + 2*f2 + 2*f3 + f4);
 return 0;
}

int main(void)
{
 int nsteps = 5000;
 float xv[nsteps];
 float yv[nsteps];
 float tv[nsteps];
 double rhoc;
 printf("Central density?");
 scanf("%lf",&rhoc);
 double coord[3] = {rhoc,0,1E-6}; /* Initial coordinates, r cannot be zero */
 xv[0] = coord[0];
 yv[0] = coord[1];
 tv[0] = coord[2];
 int i;
 double * p = &coord[0];    /* Point to first coordinate */
 for (i=1;i<nsteps;i++)
 {
  rk_coupled(p, 0.001,f,g);   /* Updates value       */
  xv[i] = coord[0];
  yv[i] = coord[1];
  tv[i] = coord[2];
  if (yv[i]<=yv[i-1])
  {
   int IER = cpgbeg(0,"/PNG",1,1);
   if (IER != 1)
   {printf("failure");
   exit(0);}
   if (yv[i-1]>rhoc)
   {
    cpgenv(0,tv[i-1]+0.5,0.0,yv[i-1],0,1);  /* If mass is higher than the central density, set the y-axis to the mass value */
   }
   else
   {
    cpgenv(0,tv[i-1]+0.5,0.0,rhoc,0,1); /* Scale the x-axis slight ab0ve the max value so the graph is clearer */
   }
   cpglab("x","y","Plot of Interior Mass (Red) and Density (White)");
   cpgline(i-1,tv,xv);     /* Plot of density. Take i-1 as i's element is the wrong one */
   cpgsci(2);
   cpgline(i-1,tv,yv);     /* Plot of interior mass, red in colour */
   cpgend();
   printf("Star's mass: %f Radius: %f\n",yv[i-1],tv[i-1]);
   exit(0);
  }
 }
}
