//
// (Extended) Mean-Field Network Simulations
//
// Lausanne, June 3rd 2008 - Michele Giugliano, PhD.
// mgiugliano@gmail.com
//
// Mac users   : download "Xcode" from Apple Web Site
// Win users   : download Cygwin from www.cygwin.com - be sure to include 'gcc compiler'
// Linux users : download distribution packages containing the 'gcc compiler'
//
// Compile with:    gcc -o meanfield source/meanfield.c -lm -O
// Matlab plot:     x = load('data.x'); figure(1); clf; subplot(2,1,1); plot(x(:,1), x(:,2)); subplot(2,1,2); plot(x(:,1), x(:,3));
//

#include <stdio.h>                      // Standard i/o library
#include <stdlib.h>                     // Standard library
#include <time.h>                       // Standard time library (for initialization purpouses)
#include <math.h>						// Standard math library
#include "../libs/rando_pyl.h"             // Needed to use gauss()    --> random number generation
#include "../libs/files_pyl.h"             // Needed by load_par()     --> i/o with parameter file
#include "../libs/nrutil.h"                // Numerical Recipes routines [custom library].
#include "../libs/nr.h"                    // Numerical Recipes routines [custom library].

#include "../libs/rando_pyl.c"             // Needed to use gauss()    --> random number generation
#include "../libs/files_pyl.c"             // Needed by load_par()     --> i/o with parameter file
#include "../libs/nrutil.c"                // Numerical Recipes routines [custom library].
#include "../libs/nr.c"                    // Numerical Recipes routines [custom library].

#define USAGE   "USAGE: %s T N C I mext sext Use\n\n"
#define MAX(a,b)    (((a)>(b)) ? (a) : (b))	// Useful macro for getting the maximum of two numbers.

#define PARFILENAME "config_files/ifparbest.par"     // Parameter file, for single-neuron properties (from Giugliano et al., 2004).

//-----------------
FILE *fopen(), *fp, *fq;				// Output file pointers.
double p[10];                           // Array containing the model parameters (5 or 6). 
double mi, si, taui;                    // Actual mean, stdev [pA] and corr. length of I [s].
double mext, sext;                      // Background activity [pA].
double t, dt, sdt, T;                   // Actual time, integration step and total simulation time [ms].
double N, C, I;                         // Network-related parameters..
double noise, R;                        // Additional state variables.
double taux, alpha, X;                  // Time constant and scaling coefficients - spike-freq adaptation.
double tauD, tauF, U;                   // Parameter for short-term depression and facilitation.
double r, u;                            // State variables for short-term depression and facilitation.
//-----------------
 
 
//----- FUNCTION PROTOTYPES ------------------------------------------------------------------------
double phi(double);                            // nerf() function, used in if_tf().
double if_tf(double *, double, double, double);// It computes the mean firing rate 
                                               // of the Leaky IF neuron without adaptation.
void load_par(char *);                         // It loads the model parameters from file.
void init();								   // Initialization of the simulation.
void print();								   // Data output routine
//---------------------------------------------------------------------------------------------------


int main(int argc, char **argv)  {				// main
 double m, s;
 int    bool  = 0;
 double tlast = -999;							// it is initialized to (almost) -infinity
 
 if (argc < 8) {								// Should the software be called with a wrong input arguments number
  printf(USAGE, argv[0]);						// information on its usage are printed on the standard output.
  exit(0);										// However, in this case the program exits.
 }
 
 init();										// The simulation is being initialized.
 
 load_par(PARFILENAME);                          // Best fitting parameters are read from file (see Giugliano et al., 2004).  
 T    = atof(argv[1]);                           // Total simulation lifetime.. [ms]. 
 N    = atof(argv[2]);                           // Size of the simulated network.
 C    = atof(argv[3]);                           // Probability of a pair-connection.
 I    = atof(argv[4]);                           // Mean of the synaptic efficacy.
 mext = atof(argv[5]);                           // Background synaptic activity [pA].
 sext = atof(argv[6]);                           // Background synaptic activity [pA].
  
 U    = atof(argv[7]); 
  
 fp = fopen("simulation_results/data.x", "w");						 // Output file is opened here..
 fq = fopen("simulation_results/bursts.x", "w");					 // Output file is opened here..
 
 printf("Mean Field Simulation - (c) 2008 Michele Giugliano, PhD.\n\n");
 printf("N = %f, C = %f, I = %f, mext = %f, sext = %f\n\n", N, C, I, mext, sext);
 
 // e.g. ./newmeanfield.exe 10000 100 .38 40 20 90 0.5
 
 //         A [pA]  U    F [ms]  D [ms]
 // control 173   0.52  26   419   
 // cnt     64   0.55   104   255   
 
 tauF = 1.;  								// Facilitating time constant [ms]
 tauD = 255.;								// Depressing time constant [ms]
 //U    = 0.55;								// Usage effective parameter

 taui = 10.;                                // Correlation time length [ms]..
 taux = 700.;                               // Spike-frequency adaptation [ms].. 

 R    = 0.;									// Mean firing rate [kHz]
 X    = 0.;									// Spike-frequency adaptation state variable..
 r    = 0.;									//
 u    = U;									//
 
 m    = mext;								//
 s    = sext*sext;							//
 
 
 while (t <= T) {   // Main simulation cycle..
  alpha =0;
//  X     +=  (R - X) * dt/taux;

  r     +=  (1. - r) * dt/tauD - u * r * R * dt;
  r     =   (r > 0) ? r : 0.;
 
  u     +=  (U  - u) * dt/tauF + U * (1. - u) * R * dt;
  u     =   (u > 0) ? u : 0.;
  u     =   (u > 1) ? 1 : u;
    
  m     += (N * C * (I*u*r/U) * R * taui - m) * dt/taui;
  s     += (0.5 * N * C * ((I*u*r/U) * (I*u*r/U)) * R * taui - s) * dt/(taui/2.);
	  
  R    += dt/(2.) * (0.001 * if_tf(p, m - alpha * X + mext, sqrt(s + sext*sext), taui) - R);
  R    = R + 1 * gauss() * sqrt(R / N);
  R    = (R > 0) ? R : 0.;
  

  if ((R > 20.*0.001) & (bool==0) & (t-tlast)>100.) {
   fprintf(fq, "%f\n", t);
   bool = 1;
   tlast=t;
  }
  if ((R < 20.*0.001) & (bool==1)) {
   fprintf(fq, "%f 0\n", t);
   bool = 0;
  }
   
  t += dt;
  
  if (fmod(t,1*dt)<=dt) { print();  fflush(NULL); }
  
 } // end while()


fclose(fp);
fclose(fq);

return 0;
} // end main()
//------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------
void init() {
 time_t t1;

 t  = 0.;                // Current time.. [ms].
 dt = 1.;                // Integration time step.. [ms].
 sdt= sqrt(dt/1000.);    // Square root of the 'dt' [ms^0.5].
  
 //alpha = 6.232447;
 (void) time(&t1); 
 srand49((long) t1);
 //printf(">>> %d\n", (long) t1);    
return;
} // end init()
//------------------------------------------------------------------------------------------------




//------------------------------------------------------------------------------------------------
void print() {
 fprintf(fp, "%f %f %f %f\n", t, 1000*R, X, r*u/U);
return;
}
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
double phi(double x) {                           // nerf() function, used in evaluating if_tf()
  return (1.772453851*nerf(x));                  // note: sqrt(pi) = 1.77245..
} // end phi()
//
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
// if_tf() computes the mean firing rate of the IF neuron without adaptation. [return Hz]
double if_tf(double *p, double mi, double si, double taui)   {
  // mu     : expected value of the experimentally injected gauss-distributed current (Ornstein-Uhlenbeck process).
  // sigma  : standard deviation of the exp. injected gauss-distributed current (Ornstein-Uhlenbeck process).
  // Tm [ms]: parameter of the IF model (it represents the passive membrane time constant).
  // Tarp [ms]: parameter of the IF model (it represents the absolute refractory period of the spike emission process).
  // theta [mV]   : parameter of the IF model (it represents the excitability threshold). (NOTE: this is fixed @ 20 mV).
  // double H [mV]: parameter of the IF model (it represents the reset hyperpolarizing voltage).
  //
  // freq [Hz]  : output mean firing rate, function of the input (mu,sigma) and the model parameters.

  // Tm    = p[1]; (if model == 1, i.e. IF has been chosen) 
  // beta  = p[1]; (if model == 0, i.e. LIF has been chosen)
  // C     = p[2];
  // H     = p[3];
  // Tarp  = p[4];
  // theta = 20.;
  // mi     = DATA[i][0]; si = DATA[i][1]; taui = DATA[i][2];

  double freq, a, b, integ, tmp, Tm, C, mu, sigma, Tarp, theta, H;

  Tm    = p[1];
  C     = p[2];
  mu    = (mi)/C;
  sigma = si/C; 
  //sigma = si*sqrt(2*taui)/C;
  H     = p[3];
  Tarp  = p[4];
  theta = 20.;

    if (sigma <= 0.) {  // if the stdev is < 5pA, don't use Ricciardi's TF, but assume sigma = 0.
     freq = 1000. / (Tarp + Tm * log((H-mu*Tm)/(theta-mu*Tm)) );
    return freq;
   }

  tmp = (sigma*sqrt(Tm));    
  a   = (H-mu*Tm)/tmp;       
  b   = (theta-mu*Tm)/tmp; // integration boundaries definition.

  if (a > 4.9 || b > 4.9)  // if the integral boundaries are too large, the int. --> +infty
   return 0.;  // freq = 0.;
  
  integ = qsimp(phi, a, b, 1.e-6);
  
  if(integ == -1.) {
   fprintf(stderr,"(if_tf(): (WARNING): 'qsimp' returned -1.. [%f %f; %f %f %f %f]\n",mi,si,Tm,C,H,Tarp); 
   return 0.; // freq = 0.; [AS A DEFAULT IN SUCH A CASE]
  }
  
  freq = 1000./(Tarp + Tm*integ);   // [Hz]
  if (freq<0.)     return 0.;
  else             return freq;
} // end if_tf()
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void load_par(char *filename)    {      // Function to read best pars from file..
  readline(filename, p);                // Par parsing [see 'files.c' for more details..].
  alpha = p[0];
  return;
} // end load_par()
//------------------------------------------------------------------------------------------------
