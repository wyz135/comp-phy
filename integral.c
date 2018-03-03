#include <stdio.h>
#include <math.h>

/* log(1+x)/x function */
double log1(double x)
{
  if (x == 0)
   {
    return 1.0;
   }
  else
   {
    return log(1+x)/x;
   }
}

/* Executes rectangular integral rule */
double rect_int(double low, double high, int nsteps, double (* func)(double))
{
    int i=0;
    double sum = 0.0;
    double h = (high - low)/((double) nsteps);
    double x= low;
    sum = func(x);
    for(i=1;i < nsteps; i++)
    {
         x = x + h;
         sum = sum + func(x);
    }
    return sum * h;
}

/* Executes the trapezoid rule on func */
double trap_int(double low, double high, int nsteps, double (* func)(double))
{
 int i=0;
 double sum = func(low) + func(high);;
 double h = (high - low)/((double) nsteps);
 double x = low;
 for(i=1; i < nsteps; i++)
 {
  x = x + h;
  sum = sum + 2* func(x);
 }
 return (sum * h)/2;
}

/* Function to execute Simpson's rule */
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

int main(void)
{
  int nsteps =1000;
  printf("Input number of steps:");
  scanf("%d",&nsteps); /* Prompt user to input number of steps */
  double pi=4.0*atan(1.0);
  double integ_actual = pi*pi/12;
  double integ_rect = rect_int(0,1,nsteps,log1);
  double error_rect = 100*fabs((integ_rect - integ_actual)/integ_actual);
  double integ_trap = trap_int(0,1,nsteps,log1);
  double error_trap = 100*fabs((integ_trap - integ_actual)/integ_actual);
  double integ_simp = simp_int(0,1,nsteps,log1);
  double error_simp = 100*fabs((integ_simp - integ_actual)/integ_actual);
  printf("Actual value: %f\n",integ_actual);
  printf("Rectangular approx: %f Relative Error: %f%%\n", integ_rect, error_rect);
  printf("Trapezoid approx: %f Relative Error: %f%%\n", integ_trap, error_trap);
  printf("Simpson approx: %f Relative Error: %f%%\n", integ_simp, error_simp);
}

