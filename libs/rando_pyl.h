//
// rando_pyl.h
// (it assumes include math.h)
//
// Available routines (from Numerical Recipes in C - Press et al. 1992):
//
// double drand49():        returns a uniform "double" random number bewteen 0.0 and 1.0
// double srand49(long):    initializes the random seed (i.e. "rand49_idum") and returns a double as drand49()
// double drand10():        returns "1" or "0" with probability 0.5 and 0.5, respectively.
// double srand10(long):    initializes the random seed (i.e. "iseed") and returns a double as drand49()
// double gauss():          returns a Gauss-distributed "double" random number, with zero mean and unitary var. 
// double RV(p):            returns "1" with probability "p".
//
// Written by Stefano Fusi, Maurizio Mattia et al. @Rome
// Used and Revised on 11 Oct 2003 by Michele Giugliano, PhD (info and bug reports to michele@giugliano.info)
//

static long rand49_idum    = -77531;     // Random seed for the random number generator "drand49()", "gauss()", "RV()".
static unsigned long iseed = 31277;      // Random seed for the random number generator "drand10()"

double drand49(void);
double srand49(long);
int    srand10(long);
int    drand10(void);
double gauss(void);
int    RV(double);
