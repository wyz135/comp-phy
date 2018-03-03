#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "cpgplot.h"

// The number of particles
#define num 10000000
// Mass of mesons and muons in GeV/c^2
#define mass_pion 0.1396
#define mass_kaon 0.4937
#define mass_muon 0.1057
// Life time of the mesons in 1E-8s
#define life_pion 2.608
#define life_kaon 1.237
// Speed of light in meters per 1E-8 seconds
#define c 2.9979

int main(void)
{
 double mean = 200; double sd = 10; /* Momentum's mean and standard deviation */
 double u1; double u2;
 static float pion_radpos[num]; static float kaon_radpos[num]; /* Initialize for radial position    */
 static float pion_neu_p[num]; static float kaon_neu_p[num];   /* Initialize for neutrino momentums */
 // Initializing seed for drand
 int n = (int) time(NULL);
 srand(n);
 double d = 2; /* Deliberately choose a number bigger than 1 to initialize */
 int i; n = 0; int m = 0;
 for (i=0;i<num;i++)
 {
  // Starts the Box-Muller method
  while (d>1)
  {
   double x1 = drand48(); double x2 = drand48();
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
   double neutrino_p = (mass_pion*mass_pion - mass_muon*mass_muon)/(2*mass_pion);
   double angle = 3.14159*drand48();       /* Generates the angle of decayed neutrino from beam axis                        */
   double pl_rest = neutrino_p*cos(angle); /* Resolves into longitudinal component in rest frame                            */
   double pt_lab = neutrino_p*sin(angle);  /* Resolves into transverse component in rest frame, which is equal to lab frame */
   // Evaluates the relativistic beta and Lorentz factor
   double beta_pion = fabs(pion_p)/sqrt(pion_p*pion_p + mass_pion*mass_pion);
   double gamma_pion = 1/sqrt(1 - beta_pion*beta_pion);
   // Transforming into lab frame and finding radial positions
   double pl_lab = gamma_pion*(pl_rest + neutrino_p*beta_pion);
   pion_neu_p[n] = pow((pl_lab*pl_lab + pt_lab*pt_lab),0.5);
   pion_radpos[n] = (700 - pion_s)*(pt_lab/pl_lab); /* Evalute the radial position */
   n++;
  }
  // After decaying, for kaons
  double k = drand48(); /* To simulate the 64% probability of Kaon decaying into neutrinos*/
  if (kaon_s < 300 & k < 0.64)
  {
   double neutrino_p = (mass_kaon*mass_kaon - mass_muon*mass_muon)/(2*mass_kaon); /* Momentum of neutrino in rest frame     */
   double angle = 3.14159*drand48();       /* Generates the angle of decayed neutrino from beam axis                        */
   double pl_rest = neutrino_p*cos(angle); /* Resolves into longitudinal component in rest frame                            */
   double pt_lab = neutrino_p*sin(angle);  /* Resolves into transverse component in rest frame, which is equal to lab frame */
   // Evaluates the relativistic beta and Lorentz factor
   double beta_kaon = fabs(kaon_p)/sqrt(kaon_p*kaon_p + mass_kaon*mass_kaon);
   double gamma_kaon = 1/sqrt(1 - beta_kaon*beta_kaon);
   // Transforming into lab frame and finding radial positions
   double pl_lab = gamma_kaon*(pl_rest + neutrino_p*beta_kaon);
   kaon_neu_p[m] = pow((pl_lab*pl_lab + pt_lab*pt_lab),0.5);
   kaon_radpos[m] = (700 - kaon_s)*(pt_lab/pl_lab); /* Evaluate the radial position */
   m++;
  }
   d = 2; /* Reset d to 2 so Box-Muller method can start again */
 }
 // Begin plotting
 printf("Number of pion's neutrino detected: %d\nNumber of Kaon's neutrino detected: %d\n", n, m*7/43);
 int IER = cpgbeg(0,"/PNG",1,1);
 if (IER != 1)
  {printf("failure");
   exit(0);}
 cpgenv(0,1.5,0,250,0,1);
 cpgpt(n,pion_radpos,pion_neu_p,1);      /* Plot for pion's neutrinos */
 cpgpt(m/6.14,kaon_radpos,kaon_neu_p,1); /* Number of points plot reduced by a factor of 6.14 to reflect the correct ratio of beam */
 cpglab("Radial Position (m)","Momentum (GeV/c)","Radial Position vs Momentum");
 cpgend();
}
