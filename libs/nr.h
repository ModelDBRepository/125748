// nr.h    ---   Oct 2000
//
// 1. nr.c is a file containing routines from NR
// 2. include nr.h in any project you want to use those routines
// 3. to add a routine it's also necessary to add the prototype in nr.h
//
// 4. if you want to use it without a makefile you have to include:
/*
  #include "nrutil.h"
  #include "nrutil.c"
  #include "nr.c"
*/
// in your program.

double nerf(double z);

double trapzd(double (*func)(double), double a, double b, int n);

double qsimp(double (*func)(double), double a, double b, double eps);

void sort(unsigned long n, float arr[]);

float probks(float alam);

// ksone: modified (23 Oct 2000) by S Fusi & G La Camera
// to allow for correlations via the parameter 'float corrn':
// corrn: number of correlated points (~tau_corr/dt) 
void ksone(float data[], unsigned long n, float corrn, float (*func)(float), float *d, float *prob);

