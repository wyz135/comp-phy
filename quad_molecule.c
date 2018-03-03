#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* Experimental parameters and global variables */

double xpi = 3.14159265359;
double V0 = 4.747;             /* Measured value of V0 in eV              */
double req = 0.74166;          /* Measured equilibrium length in angstrom */          
double EnergyLevels[15] = {-4.477, -3.962, -3.475, -3.017, -2.587, -2.185, -1.811, -1.466, -1.151, -0.867, -0.615, -0.400, -0.225, -0.094, -0.017};                       /* Experimentally determined energy levels */
double vgamma = 31.0556;       /* Gamma constant in adjusted units        */

double simp_int(double low, double high, int nsteps, double (* func)(double))
{
 int i=0;
 double sum = func(low) + func(high);
 double h = (high - low)/((double) nsteps);
 double x = low;
 for(i=1; i < nsteps; i++)
 {
  x = x + h;
  if((i/2)*2!=i) /* If i is odd */
   {sum = sum + 4*func(x);}
  else           /* If i is even */
   {sum = sum + 2*func(x);}
  }
 return (sum * h)/3;
}

double find_root(double guess, double prec, double fstep, double (* func)(double))
{
 int i = 0;
 int max_step = 10000;
 double xprev,xn,xnext;
 double yprev,yn,ynext;
 xprev = guess;
 yprev = func(xprev);
 xnext = xprev + fstep;
 ynext = func(xnext);
 /*Start Scanning so yprev and ynext has opposite signs */
 for(i=0;yprev*ynext > 0 && i < max_step;i++) /* i < 10000 To prevent infinite loop */
  {xprev = xnext;
   yprev = ynext;
   xnext = xnext + fstep;
   ynext = func(xnext);
  }
 if (i == max_step) /* Terminates the program if there is no crossover */
 {
  printf("Fail to find zero between %f and %f, algorithm terminated.\n", guess, guess + max_step*fstep);
  exit(0);
 }
 /* Now start the False-position algorithm */
 i = 0; /* Reset i */
 yn = prec + 1;
 for(i=0;fabs(yn)>prec && i < max_step;i++)
  {
   xn = xnext - (ynext*(xnext - xprev)/(ynext - yprev));
   yn = func(xn);
   if (yprev*yn < 0)
   {
    xnext = xn;
    ynext = yn;
   }
   else
   {
    xprev = xn;
    yprev = yn;
   }
  }
 return xn;
}

/* Algorithm to evaluate the position of the lowest element in the array */
int minimum(double array[1000])
{
 int i=1;
 int min = 0;
 double x = array[0];
 while (i<1000)
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

/* Begin the program to first optimize the value of a, then find the energy levels */
int main(void)
{
 double energy[15];
 double Error[1000];
 double a = 0.2; /* Initialize the value of a */
 
 double quad_integ(double en, double x) /* The integrand of (1) for quadratic potential */
 {
  return sqrt(fabs((en - 4*((x/a)-1)*((x/a)-2)))); /* fabs in case it is negative */
 }
 
 double quad_action(double en)
 {
  double rmin = (a*(3 - sqrt(1+en)))/2;
  double rmax = (a*(3 + sqrt(1+en)))/2;
  double quad_integfixed(double x) /* Specifying the value of en */
  {
   return quad_integ(en, x);
  }
  return simp_int(rmin, rmax, 10000, quad_integfixed);
 }
 
 /* Equations to evaluate the first 2 energy levels */
 double rooteq0(double en)
 {
  return quad_action(en) - (0.5)*xpi/(vgamma*sqrt(V0));
 }
 double rooteq1(double en)
 {
  return quad_action(en) - (1.5)*xpi/(vgamma*sqrt(V0));
 }

 /* Evaluate the errors for different values of a */
 int n;
 for (n=0;n<1000;n++)
 {
  double energy0 = find_root(-1.0, 1E-6, 0.1, rooteq0);
  double energy1 = find_root(-1.0, 1E-6, 0.1, rooteq1);
  double error = pow((energy0*V0 - EnergyLevels[0]),2) + pow((energy1*V0 - EnergyLevels[1]),2);
  Error[n] = error;
  a = a + 0.001;
 }
 int min = minimum(Error); /* Determines which member is the lowest */
 a = a - 0.001*(1000-min); /* Adjusts back the value of a           */
 printf("Optimal value of a: %f\n", a);
 /* Print out the evaluated values of energy levels given using the optimzed a */
 printf("Quadratic\n");
 printf("Numerical  Analytical\n");
 for (n=0; n<15; n++)
 {
  double rooteq(double en)
  {
   return quad_action(en) - (n+0.5)*xpi/(vgamma*sqrt(V0));
  }
  energy[n] = find_root(-1.0, 0.0000001, 0.01, rooteq);
  double analytic_energy = -1 + 4*(n+0.5)/(vgamma*sqrt(V0)*a);
  printf("%f,%f\n", energy[n]*V0, analytic_energy*V0);
 }
}
