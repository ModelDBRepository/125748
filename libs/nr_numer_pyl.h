//
// nr_numer.h
//
// Written by Stefano Fusi, Maurizio Mattia et al. @Rome
// Used and Revised on 13 Oct 2003 by Michele Giugliano, PhD (info and bug reports to michele@giugliano.info)
//

float *vector();
float **matrix();
float **convert_matrix();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
double **submatrix();
void free_vector();
void free_dvector();
void free_ivector();
void free_matrix();
void free_dmatrix();
void free_imatrix();
void free_submatrix();
void free_convert_matrix();
void nrerror();
void ludcmp(double **a,int n,int *indx,double *d);
void lubksb(double **a,int n,int *indx,double b[]);
double nerf(double z);
double trapzd(double (*func)(double),double a,double b,int n);
double qsimp(double (*func)(double),double a,double b);
void balanc(double **a,int n);
void elmhes(double **a,int n);
void hqr(double **a,int n,double wr[],double wi[]);
int Active_bits(int i, int p);
float gammln(float xx);
float factln(int n);
float bico(int n,int k);
