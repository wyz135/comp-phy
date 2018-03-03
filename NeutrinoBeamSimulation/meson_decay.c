#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "cpgplot.h"

// Number of particles
#define num 10000000
// Meson's mass in GeV/c^2
#define mass_pion 0.1396
#define mass_kaon 0.4937
// Meson's lifetime in 1E-8s
#define life_pion = 2.608
#define life_kaon = 1.237
// Speed of light in meter per 1E-8 seconds
#define c = 2.9979

int main(void)
{
 double mean = 200; double sd = 10; /* Momentum's mean and standard deviation */
 double u1; double u2;
 /* Initializing seed for drand */
 int n = (int) time(NULL);
 srand(n);
 static float pion[num];
 static float kaon[num];
 double d = 2; /* Deliberately choose a number bigger than 1 to initialize */
 int i;
 for (i=0;i<num;i++)
 {
  /* Starts the Box-Muller method to generate the meson's momentum */
  while (d>1)
  {
   double x1 = drand48();
   double x2 = drand48();
   u1 = 2*x1 -1;
   u2 = 2*x2 -1;
   d = u1*u1 + u2*u2;
  }
  double ln = -log(d);
  double p1 = sd*u1*sqrt(ln/d) + mean;
  double p2 = sd*u2*sqrt(ln/d) + mean;
  // p1 and p2 allocated as pion's and kaon's momentum respectively
  pion[i]=p1; kaon[i]=p2;
  d = 2;
 }
 int IER = cpgbeg(0,"/PNG",1,1);
 if (IER != 1)
  {printf("failure");
   exit(0);}
 int bin = 200; /* Bin cannot exceed 200 */
 // Plotting histograms for the mesons' momentum
 cpghist(num,pion,mean-4*sd,mean+4*sd,bin,0);
 cpglab("Momentum (GeV/c)","Count","Pion's Momentum");
 cpgpage();
 cpghist(num,kaon,mean-4*sd,mean+4*sd,bin,0);
 cpglab("Momentum (GeV/c)","Count","Kaon's Momentum");
 cpgpage();
 for (i=0; i<num; i++)
 {
  double t1 = drand48(); double t2 = drand48(); /* Generates another random number for decay distances */
  pion[i] = - ((pion[i]*life_pion*c/mass_pion)*log(t1));
  kaon[i] = - ((kaon[i]*life_kaon*c/mass_kaon)*log(t2));
 }
 // Plotting histograms for the mesons' decay distances
 cpghist(num,pion,0,30000,bin,0);
 cpglab("Distance (m)","Count","Pion's decay");
 cpghist(num,kaon,0,20000,bin,0);
 cpglab("Distance (m)","Count","Kaon's decay");
 cpgend();
}
