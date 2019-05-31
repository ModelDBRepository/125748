//
// rando_pyl.c
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



//---------------------------------------------------------------
#define M 714025
#define IA 1366
#define IC 150889
double drand49()
{
  static long iy, ir[98];
  static int iff = 0;
  int j;
  
  if (rand49_idum < 0 || iff == 0)  {
    iff = 1;
    if ((rand49_idum = (IC - rand49_idum) % M) < 0)
      rand49_idum = (-rand49_idum);
    for (j = 1; j <= 97; j++)      {
      rand49_idum = (IA * (rand49_idum) + IC) % M;
      ir[j] = (rand49_idum);
    }
    rand49_idum = (IA * (rand49_idum) + IC) % M;
    iy = (rand49_idum);
  }
  j = 1 + 97.0 * iy / M;
  if (j > 97 || j < 1)
    printf ("RAN2: This cannot happen.");
  iy = ir[j];
  rand49_idum = (IA * (rand49_idum) + IC) % M;
  ir[j] = (rand49_idum);
  return (double) iy / M;
} // end drand49()
//---------------------------------------------------------------
double srand49(long seme)
{
  rand49_idum = (-seme);
  return drand49 ();
} // end srand49()

#undef M
#undef IA
#undef IC
//---------------------------------------------------------------



//---------------------------------------------------------------
#define IB1 1
#define IB2 2
#define IB5 16
#define IB18 131072
#define MASK IB1+IB2+IB5

int srand10(long seme)
{
  iseed = seme;
} // end srand10()
int drand10()
{
  if (iseed & IB18)    {
    iseed = ((iseed ^ MASK) << 1) | IB1;
    return 1;
  }
  else     {
    iseed <<= 1;
    return 0;
  }
} // end drand10
#undef MASK
#undef IB18
#undef IB5
#undef IB2
#undef IB1
//---------------------------------------------------------------



//---------------------------------------------------------------
double gauss()
{
  static int iset = 0;
  static double gset;
  double fac, r, v1, v2;
  
  if (iset == 0)    {
    do   {
      v1 = 2.0 * drand49 () - 1.0;
      v2 = 2.0 * drand49 () - 1.0;
      r = v1 * v1 + v2 * v2;
    }
    while (r >= 1.0);
    fac = sqrt(-2.0 * log (r) / r);
    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  }
  else
    {
      iset = 0;
      return gset;
    }
} // end gauss()
//---------------------------------------------------------------



//---------------------------------------------------------------
int RV(double p)
{
  return (drand49()<p);
} // end RV()
//---------------------------------------------------------------
