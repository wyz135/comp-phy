#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "cpgplot.h"

// Number of particles
#define num 1000000
// Meson's mass in GeV/c^2
#define mass_pion 0.1396
#define mass_kaon 0.4937
#define mass_muon 0.1057
// Meson's lifetime in 1E-8s
#define life_pion 2.608
#define life_kaon 1.237
// Speed of light in meter per 1E-8 seconds
#define c 2.9979

int main(void)
{
 double mean = 200; double sd = 10;                            /* Momentum's mean and standard deviation */
 double u1; double u2;
 static float pion_neu_p[num]; static float kaon_neu_p[num];   /* Initialize arrays for pion and kaon momentum in lab frame */
 /* Initializing seed for drand */
 int n = (int) time(NULL);
 srand(n);
 double d = 2; /* Deliberately choose a number bigger than 1 to begin Box-Muller method */
 int i; n = 0; int m = 0;
 for (i=0;i<num;i++)
 {
  /* Starts the Box-Muller method */
  while (d>1)
  {
   double x1 = drand48(); double x2 = drand48(); /* Generates two random numbers for meson's momentum */
   u1 = 2*x1 -1;
   u2 = 2*x2 -1;
   d = u1*u1 + u2*u2;
  }
  double ln = -log(d);
  double pion_p = sd*u1*sqrt(ln/d) + mean;                    /* Evaluates the momentum of pion                      */
  double kaon_p = sd*u2*sqrt(ln/d) + mean;                    /* Evaluates the momentum of kaon                      */
  double t1 = drand48(); double t2 = drand48();               /* Generates another random number for decay distances */
  double pion_s = - ((pion_p*life_pion*c/mass_pion)*log(t1)); /* Evaluates the decay distance of pion                */
  double kaon_s = - ((kaon_p*life_kaon*c/mass_kaon)*log(t2)); /* Evaluates the decay distance of kaon                */
  // After decaying, only neutrino produced before 300m will be detected
  // After decaying, for pions
  if (pion_s < 300)
  {
   double neutrino_p = (mass_pion*mass_pion - mass_muon*mass_muon)/(2*mass_pion); /* Momentum of neutrino in rest frame     */
   double angle = 3.14159*drand48();       /* Generates the angle of decayed neutrino                                       */
   double pl_rest = neutrino_p*cos(angle); /* Resolves into longitudinal component in rest frame                            */
   double pt_lab = neutrino_p*sin(angle);  /* Resolves into transverse component in rest frame, which is equal to lab frame */
   // Evaluates the relativistic beta and Lorentz factor
   double beta_pion = fabs(pion_p)/sqrt(pion_p*pion_p + mass_pion*mass_pion);
   double gamma_pion = 1/sqrt(1 - beta_pion*beta_pion);
   // Transforming into lab frame and finding neutrino's momentum
   double pl_lab = gamma_pion*(pl_rest + neutrino_p*beta_pion);
   pion_neu_p[n] = pow((pl_lab*pl_lab + pt_lab*pt_lab),0.5);
   n++;
  }
  //  After decaying, for kaons
  double k = drand48(); /* To simulate the 64% probability of Kaon decaying into neutrinos */
  if (kaon_s < 300 & k < 0.64)
  {
   double neutrino_p = (mass_kaon*mass_kaon - mass_muon*mass_muon)/(2*mass_kaon); /* Momentum of neutrino in rest frame     */
   double angle = 3.14159*drand48();       /* Generates the angle of decayed neutrino                                       */
   double pl_rest = neutrino_p*cos(angle); /* Resolves into longitudinal component in rest frame                            */
   double pt_lab = neutrino_p*sin(angle);  /* Resolves into transverse component in rest frame, which is equal to lab frame */
   // Evaluates the relativistic beta and Lorentz factor
   double beta_kaon = fabs(kaon_p)/sqrt(kaon_p*kaon_p + mass_kaon*mass_kaon);
   double gamma_kaon = 1/sqrt(1 - beta_kaon*beta_kaon);
   // Transforming into lab frame and finding neutrino's momentum
   double pl_lab = gamma_kaon*(pl_rest + neutrino_p*beta_kaon);
   kaon_neu_p[m] = pow((pl_lab*pl_lab + pt_lab*pt_lab),0.5);
   m++;
  }
   d = 2; /* Reset d to 2 so Box-Muller method can start again */
 }
 // Begin plotting
 int IER = cpgbeg(0,"/PNG",1,1);
 if (IER != 1)
  {printf("failure");
   exit(0);}
 cpghist(n,pion_neu_p,0,100,200,0); /* Plot for pion's neutrinos */
 cpglab("Momentum (GeV/c)", "Count", "Pion's momentum in lab frame");
 cpghist(m/6.14,kaon_neu_p,0,250,200,0); /* Number of points plot reduced by a factor of 6.14 to reflect the correct ratio of beam */
 cpglab("Momentum (GeV/c)", "Count", "Kaon's momentum in lab frame");
 cpgend();
}
