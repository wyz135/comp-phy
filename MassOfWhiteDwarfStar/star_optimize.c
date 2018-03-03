#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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

/* Algorithm to evaluate the position of the lowest element in the array */
int minimum(double array[5000])
{
 int i=1;
 int min = 0;
 double x = array[0];
 while (i<5000)
  {
   if (array[i] < x)
   {
    x = array[i];
    min = i;
    i++;
   }
   else
   {
    i++;
   }
  }
 return min;
}

int main(void)
{
 int nsteps = 5000; /* Now used to specify array size */ 
 double mass[nsteps];
 double radius[nsteps];
 float xv[nsteps];
 float yv[nsteps];
 float tv[nsteps];
 double error[nsteps];
 double Ye;
 double actual_mass;
 double actual_radius;
 printf("Ye value?\n");
 scanf("%lf",&Ye);
 printf("Mass of star?\n");
 scanf("%lf",&actual_mass);
 printf("Radius of star?\n");
 scanf("%lf",&actual_radius);
 double rhoc[nsteps];
 rhoc[0] = 0.01;
 int n;
 for (n=1;n<nsteps;n++)
 {
  rhoc[n] = rhoc[n-1]+0.01;
 }
 for (n=0; n<nsteps; n++)
 {
  double coord[3] = {rhoc[n],0,1E-6};                         /* Initial conditions, r cannot be zero */
  xv[0] = coord[0]; yv[0] = coord[1]; tv[0] = coord[2];
  double * p = &coord[0];                                     /* Point to first coordinate */
  rk_coupled(p, 0.001,f,g);                                   /* Updates value       */
  xv[1] = coord[0]; yv[1] = coord[1]; tv[1] = coord[2];       /* Assign 2nd elements in the array for comparison in the loop */
  int i;
  for (i=2;yv[i-1]>yv[i-2] && i<nsteps ;i++)                               /* Immediately terminates loop when the highest mass is reached */
  {
   rk_coupled(p, 0.001,f,g);   /* Updates value       */
   xv[i] = coord[0];
   yv[i] = coord[1];
   tv[i] = coord[2];
  }
  mass[n] = yv[i-2]*5.67*Ye*Ye/1.989; /* Evaluates mass and converts to solar units */
  radius[n] = tv[i-2]*7.72*Ye/(695.7); /* Evaluates radius and converts to solar units */
  error[n] = fabs((mass[n] - actual_mass)/actual_mass) + fabs((radius[n] - actual_radius)/actual_radius); /* Computes relative error */
 }
 int min = minimum(error); /* Identifies the position of minimum error */
 printf("Central Density: %f Mass: %f Radius: %f Error:%f\n",rhoc[min], mass[min], radius[min],error[min]);
}
