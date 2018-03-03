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

int euler_coupled(double * coord, double step, double (* f) (double,double,double), double (* g) (double,double,double))
{
 /* Updates the coordinates by 1 step via Euler. Important note: For the first parameter,
    need to point coord to the first coordinate first */
 double xn = *coord;
 coord++;
 double yn = *coord;
 coord++;
 double tn = *coord;
 /* Finding the gradients */
 double f1 = f(xn,yn,tn);
 double g1 = g(xn,yn,tn);
 * coord = tn + step;
 coord = coord - 1;
 * coord = yn + step*g1;
 coord = coord - 1;
 * coord = xn + step*f1;
 return 0;
}

int rk_coupled(double * coord, double step, double (* f) (double,double,double), double (* g) (double,double,double))
{
 /* Updates the coordinates by 1 step via Runge-Kutta. Important note: FOr the first parameter,
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
 int nsteps = 500;
 float xv[nsteps];
 float yv[nsteps];
 float tv[nsteps];
 double rhoc = 10.0;              /* Central density                       */
 double h = 0.01;
 double coord[3] = {rhoc,0,1E-6}; /* Initial coordinates, r cannot be zero */
 xv[0] = coord[0];
 yv[0] = coord[1];
 tv[0] = coord[2];
 int i;
 double * p = &coord[0];    /* Point to first coordinate */
 for (i=1;i<nsteps;i++)
 {
  euler_coupled(p,h,f,g);   /* Updates value       */
  xv[i] = coord[0];
  yv[i] = coord[1];
  tv[i] = coord[2];
 }
 int IER = cpgbeg(0,"/PNG",1,1);
 if (IER != 1)
  {printf("failure");
  exit(0);}
 cpgenv(0,5.0,0.0,10.0,0,1);
 cpglab("x","y","Plot of Interior Mass and Density using Euler (red) and Runge-Kutta (green)");
 cpgsci(2);              /* Red for Euler         */
 cpgline(i-1,tv,xv);     /* Plot of density       */
 cpgline(i-1,tv,yv);     /* Plot of interior mass */
 printf("Euler's mass: %f", coord[1]);
 /* Reset the initial coordinates to use Runge-Kutta method */
 coord[0] = rhoc; coord[1] = 0; coord[2] = 1E-6;
 xv[0] = coord[0];
 yv[0] = coord[1];
 tv[0] = coord[2];
 for (i=1;i<nsteps;i++)
 {
  rk_coupled(p,h,f,g);   /* Updates value       */
  xv[i] = coord[0];
  yv[i] = coord[1];
  tv[i] = coord[2];
 }
 printf("Runge-Kutta's mass: %f\n", coord[1]);
 cpgsci(3);              /* Green for Runge-Kutta */
 cpgline(i-1,tv,xv);     /* Plot of density       */
 cpgline(i-1,tv,yv);     /* Plot of interior mass */
 cpgend();
}
