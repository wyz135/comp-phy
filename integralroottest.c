#include <stdio.h>
#include <math.h>

double find_root(double guess, double prec, double fstep, double (* func)(double))
{
 int i = 0;
 double xprev,xn,xnext;
 double yprev,yn,ynext;
 xprev = guess;
 yprev = func(xprev);
 xnext = xprev + fstep;
 ynext = func(xnext);
 /*Start Scanning so yprev and ynext has opposite signs */
 for(i=0;yprev*ynext > 0 && i < 10000;i++) /* i < 10000 To prevent infinite loop */
  {xprev = xnext;
   yprev = ynext;
   xnext = xnext + fstep;
   ynext = func(xnext);
  }
 /* Now start the False-position algorithm */
 xn = xnext - (ynext*(xnext - xprev)/(ynext - yprev));
 yn = func(xn);
 for(i=0;fabs(yn)>prec && i < 10000;i++)
  {xprev = xnext;
   yprev = ynext;
   xnext = xn;
   ynext = yn;
   xn = xnext - (ynext*(xnext - xprev)/(ynext - yprev));
   yn = func(xn);
  }
 return xn;
}

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

/* Function of x squared to be integrated */
double xsq(double x)
{
 return x*x;
}

int main(void)
{
 double guess;
 int nsteps;
 double prec;
 printf("Guess lower bound of root?\n");
 scanf("%lf", &guess);
 printf("Number of steps?\n");
 scanf("%d",&nsteps);
 printf("Root precision?\n");
 scanf("%lf",&prec);
 double integ(double x)
 {
  return simp_int(0,x,nsteps,xsq) - x;
 }
 double root = find_root(guess,prec,1.0,integ);
 printf("Root: %f\n", root);
}
