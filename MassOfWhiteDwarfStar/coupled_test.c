#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cpgplot.h"

double f(double x, double y, double t)
{
 /* DE for x */
 return -y;
}

double g(double x, double y, double t)
{
 /* DE for y */
 return x;
}

int euler_coupled(double * coord, double step, double (* f) (double,double,double), double (* g) (double,double,double))
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
 * coord = tn + step;
 coord = coord - 1;
 * coord = yn + step*g1;
 coord = coord - 1;
 * coord = xn + step*f1;
 return 0;
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
 /* Calls the pgplot function */
 int IER = cpgbeg(0,"/PNG",1,1);
 if (IER != 1)
 {
  printf("failure");
  exit(0);
 }
 double range = 10.0;
 cpgenv(0,range,-2,2,0,1);
 cpglab("x","y","Euler vs. Runge-Kutta");
 double stepsize;
 printf("Step size?");
 scanf("%lf",&stepsize);
 int nsteps = range/stepsize;
 float xv[nsteps];
 float actual_xv[nsteps];
 float yv[nsteps];
 float actual_yv[nsteps];
 float tv[nsteps];
 double coord[3] = {1.0,1.0,0.0}; /* Initial coordinates */
 xv[0] = coord[0];
 yv[0] = coord[1];
 tv[0] = coord[2];
 int i;
 double * p = &coord[0];            /* Point to first coordinate */
 for (i=1;i<nsteps;i++)
 {
  euler_coupled(p, stepsize,f,g);    /* Updates value             */
  xv[i] = coord[0];
  yv[i] = coord[1];
  tv[i] = coord[2];
 }
 for (i=0;i<nsteps;i++)
 {
  actual_xv[i] = -sin(tv[i]) + cos(tv[i]);
 }
 /* Begin plotting  */
 cpgline(nsteps,tv,xv); /* Plots for euler's method, white in colour */
 coord[0] = 1.0; coord[1] = 1.0;coord[2] = 0.0; /* Re-initializes intial coordinates */
 for (i=1;i<nsteps;i++)
 {
  rk_coupled(p, stepsize,f,g);    /* Updates value             */
  xv[i] = coord[0];
  yv[i] = coord[1];
  tv[i] = coord[2];
 }
 cpgsci(3);
 cpgline(nsteps,tv,xv); /* Plots for Runge-Kutta's method, green in colour */
 cpgsci(2);
 cpgline(nsteps,tv,actual_xv); /* Plots for actual, red in colour */
 cpgend();
}
