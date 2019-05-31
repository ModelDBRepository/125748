//
// nr_numer.c
// (it assumes include math.h)
//
// Written by Stefano Fusi, Maurizio Mattia et al. @Rome
// Used and Revised on 13 Oct 2003 by Michele Giugliano, PhD (info and bug reports to michele@giugliano.info)
//

#include <stdio.h>
#include <malloc.h>

extern int return_value;     // global error variable defined in _espo.c


void nrerror(char error_text[])
{
    fprintf(stderr,"  Numerical Recipes run-time error...\n");
    fprintf(stderr,"  %s\n",error_text);
    fprintf(stderr,"  ...now exiting to system...\n\n");
    return_value = 0;
    return;
}



float *vector(int nl, int nh)
{
    float *v;
    v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
    if (!v) nrerror("allocation failure in vector()");
    return v-nl;
}



int *ivector(int nl, int nh)
{
    int *v;

    v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
    if (!v) nrerror("allocation failure in ivector()");
    return v-nl;
}



double *dvector(int nl, int nh)
{
    double *v;

    v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
    if (!v) nrerror("allocation failure in dvector()");
    return v-nl;
}



float **matrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    float **m;

    m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m -= nrl;

    for(i=nrl;i<=nrh;i++) {
        m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
        if (!m[i]) nrerror("allocation failure 2 in matrix()");
        m[i] -= ncl;
    }
    return m;
}



double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    double **m;

    m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
    if (!m) nrerror("allocation failure 1 in dmatrix()");
    m -= nrl;

    for(i=nrl;i<=nrh;i++) {
        m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
        if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
        m[i] -= ncl;
    }
    return m;
}



int **imatrix(int nrl, int nrh, int ncl, int nch)
{
    int i,**m;

    m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
    if (!m) nrerror("allocation failure 1 in imatrix()");
    m -= nrl;

    for(i=nrl;i<=nrh;i++) {
        m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
        if (!m[i]) nrerror("allocation failure 2 in imatrix()");
        m[i] -= ncl;
    }
    return m;
}



double **submatrix(double **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl,int newcl)
{
    int i,j;
    double **m;

    m=(double **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(double*));
    if (!m) nrerror("allocation failure in submatrix()");
    m -= newrl;

    for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

    return m;
}



void free_vector(float *v, int nl, int nh)
{
    free((char*) (v+nl));
}



void free_ivector(int *v, int nl, int nh)
{
    free((char*) (v+nl));
}



void free_dvector(double *v, int nl, int nh)
{
    free((char*) (v+nl));
}



void free_matrix(float *m, int nrl, int nrh, int ncl, int nch)
{
    int i;

    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
    free((char*) (m+nrl));
}



void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
    int i;

    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
    free((char*) (m+nrl));
}



void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
    int i;

    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
    free((char*) (m+nrl));
}



void free_submatrix(b,nrl,nrh,ncl,nch)
double **b;
int nrl,nrh,ncl,nch;
{
    free((char*) (b+nrl));
}



float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
int nrl,nrh,ncl,nch;
{
    int i,j,nrow,ncol;
    float **m;

    nrow=nrh-nrl+1;
    ncol=nch-ncl+1;
    m = (float **) malloc((unsigned) (nrow)*sizeof(float*));
    if (!m) nrerror("allocation failure in convert_matrix()");
    m -= nrl;
    for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
    return m;
}



void free_convert_matrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
    free((char*) (b+nrl));
}






     /**********************************************************
      *                                                        *
      *                       ludcmp                           *
      *                                                        *
      *                                                        *
      *                                                        *
      **********************************************************/ 


#include <stdlib.h>
#include <math.h>
//#include "nrutil.h"
#include <stdio.h>


#define TINY 1.0e-20;

extern int return_value;

void ludcmp(a,n,indx,d)
int n,*indx;
double **a,*d;
{
    int i,imax,j,k;
    double big,dum,sum,temp;
    double *vv;
    void free_dvector();
    vv=dvector(1,n);
    *d=1.0;
    for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
            if ((temp=fabs(a[i][j])) > big) big=temp;
        if (big == 0.0) {
          printf(">>> Singular matrix in routine LUDCMP\n");
          return_value = 0;
          free_dvector(vv,1,n);
          return;
        }
        vv[i]=1.0/big;
    }
    for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
            sum=a[i][j];
            for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<=n;i++) {
            sum=a[i][j];
            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=1;k<=n;k++) {
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            *d = -(*d);
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        if (j != n) {
            dum=1.0/(a[j][j]);
            for (i=j+1;i<=n;i++) a[i][j] *= dum;
        }
    }
    free_dvector(vv,1,n);
}

#undef TINY




     /**********************************************************
      *                                                        *
      *                       lubksb                           *
      *                                                        *
      *                                                        *
      *                                                        *
      **********************************************************/ 


void lubksb(a,n,indx,b)
double **a,b[];
int n,*indx;
{
    int i,ii=0,ip,j;
    double sum;

    for (i=1;i<=n;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
            for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum) ii=i;
        b[i]=sum;
    }
    for (i=n;i>=1;i--) {
        sum=b[i];
        for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}





     /**********************************************************
      *                                                        *
      *                         nerf(x):                       *
      *                                                        *
      *                      x^2*(1+erf(x))                    *
      *                                                        *
      * (Qualche matematico geniale ha trovato i coefficienti  *
      *  giusti, per non dire magici...)                       *
      *                                                        *
      **********************************************************/ 


#include <math.h>
#define a1  -1.26551223
#define a2  1.00002368
#define a3  .37409196
#define a4  .09678418
#define a5  -.18628806
#define a6  .27886087
#define a7  -1.13520398
#define a8 1.48851587
#define a9 -.82215223
#define a10 .17087277

double nerf(double z)
{      
  double t,ef,at;
  double w;
  w = fabs(z);
  t = 1./(1. + 0.5 * w);
  at=a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*(a9+t*a10))))))));
  ef=t*exp(at);
  if(z>0.)
    ef = 2.*exp(w*w)-ef;
  return ef;
}




              /************************\
               *    ROUTINE trapzd    *
              \************************/

#include <stdlib.h>
#define FUNC(x) ((*func)(x))


double trapzd(double (*func)(double),double a,double b,int n)
{
  double x,tnm,sum,del;
  static double s;
  int it,j;


  if(n==1){
    return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
  } else {
    for(it=1,j=1;j<n-1;++j) it<<=1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for(sum=0.0,j=1;j<=it;++j,x+=del) sum+=FUNC(x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}


            /*********************\
             *   ROUTINE qsimp   *
            \*********************/

#include <math.h>
#include <stdio.h>
#define EPS 1.0e-6   // 1e-6
#define JMAX 20      // 20


double qsimp(double (*func)(double),double a,double b)
{
  double trapzd(double (*func)(double),double a,double b,int n);
  int j;
  double s,st,ost,os;

  //printf("\n");//
  ost=os=-1.0e30;
  for(j=1;j<=JMAX;++j){
    st=trapzd(func,a,b,j);
    s=(4.0*st-ost)/3.0;
    if(j>5)
      if(fabs(s-os)<EPS*fabs(os)||(s==0.0 && os==0.0)) return s;
    os=s;
    ost=st;
    
    //printf("--- # qsimp step: %d\r",j); //
    //fflush(stdout); //
  }

  printf("\n\nError: too many steps in routine qsimp.\n");
  return -1;   // error case
}


           /************************************************
           *                                              *
           *                     balanc                   *
           *                                              *
           * Given a matrix a[1...n][1...n] this routine  *
           * replaces it by a balanced matrix with        *
           * identical eigenvalues (by Numerical Recipes) *
           *                                              *
           ************************************************/
        

#include <math.h>
#define RADIX 2.0

void balanc(double **a,int n)
{
    int last,j,i;
    double s,r,g,f,c,sqrdx;

    sqrdx=RADIX*RADIX;
    last=0;
    while (last == 0){
      last=1;
      for (i=1;i<=n;i++) {
        r=c=0.0;
        for (j=1;j<=n;j++)
          if (j != i) {
        c += fabs(a[j][i]);
        r += fabs(a[i][j]);
          }
        if (c && r) {
          g=r/RADIX;
          f=1.0;
          s=c+r;
          while (c<g) {
        f *= RADIX;
        c *= sqrdx;
          }
          g=r*RADIX;
          while (c>g) {
        f /= RADIX;
        c /= sqrdx;
          }
          if ((c+r)/f < 0.95*s) {
        last=0;
        g=1.0/f;
        for (j=1;j<=n;j++) a[i][j] *= g;
        for (j=1;j<=n;j++) a[j][i] *= f;
          }
        }
      }
    }
}


#undef RADIX







     /******************************************************************
      *                                                                *
      *                           elmhes                               *
      *                                                                *
      * Reduction to Hessenberg form by the elimination method.        *
      * The REAL, nonsymmetric matrix a[1...n][1...n] is replaced by   *
      * an upper Hessenberg matrix with identical eigenvalues. [...]   * 
      * (By Num.Rec.)                                                  *
      ******************************************************************/




#include <math.h>

#define SWAP(g,h) {y=(g);(g)=(h);(h)=y;}

void elmhes(double **a,int n)
{

    int m,j,i;
    double y,x;

    for (m=2;m<n;m++) {
      x=0.0;
      i=m;
      for (j=m;j<=n;j++) {
        if (fabs(a[j][m-1]) > fabs(x)) {
          x=a[j][m-1];
          i=j;
        }
      }
      if (i != m) {
        for (j=m-1;j<=n;j++) SWAP(a[i][j],a[m][j])
                   for (j=1;j<=n;j++) SWAP(a[j][i],a[j][m])
                   }
      if (x) {
        for (i=m+1;i<=n;i++) {
          if (y=a[i][m-1]) {
        y /= x;
        a[i][m-1]=y;
        for (j=m;j<=n;j++)
          a[i][j] -= y*a[m][j];
        for (j=1;j<=n;j++)
          a[j][m] += y*a[j][i];
          }
        }
      }
    }
}






          /***********************************************************
           *                                                         *
           *                           hqr                           *
           *                                                         *
           * Finds all eigenvalues of an upper Hessenberg matrix     *
           * a[1...n][1...n]. On input a can be exactly as output    *
           * from elmhes(); on output it is destroyed.               *    
           * The real and imaginary parts of the eigenvalues are     *
           * returned in wr[1...n] and wi[1...n] respectively.       *
           *                      (by Num.Rec.)                      *
           ***********************************************************/              



#include <math.h>
#include <stdio.h>

#define SIGN(a,b) ((b) > 0 ? fabs(a) : -fabs(a))


extern int return_value;     // defined in espo.c


void hqr(a,n,wr,wi)
double **a,wr[],wi[];
int n;
{

        int nn,m,l,k,j,its,i,mmin;
    double s;  // float?
    double z,y,x,w,v,u,t,r,q,p,anorm;

    anorm=fabs(a[1][1]);
    for (i=2;i<=n;i++)
        for (j=(i-1);j<=n;j++)
            anorm += fabs(a[i][j]);
    nn=n;
    t=0.0;
    while (nn >= 1) {
      its=0;
      do {
        for (l=nn;l>=2;l--) {
          s=fabs(a[l-1][l-1])+fabs(a[l][l]);
          if (s == 0.0) s=anorm;
                // questo float lo lascio cosi'? NO
          if ((fabs(a[l][l-1]) + s) == s) break;
        }
        x=a[nn][nn];
        if (l == nn) {
          wr[nn]=x+t;
          wi[nn--]=0.0;
        } else {
          y=a[nn-1][nn-1];
          w=a[nn][nn-1]*a[nn-1][nn];
          if (l == (nn-1)) {
        p=0.5*(y-x);
        q=p*p+w;
        z=sqrt(fabs(q));
        x += t;
        if (q >= 0.0) {
          z=p+SIGN(z,p);
          wr[nn-1]=wr[nn]=x+z;
          if (z) wr[nn]=x-w/z;
          wi[nn-1]=wi[nn]=0.0;
        } else {
          wr[nn-1]=wr[nn]=x+p;
          wi[nn-1]= -(wi[nn]=z);
        }
        nn -= 2;
          } else {
        if (its == 160) {  // 30
          printf(">>> WARNING: Too many iterations in routine hqr(). Unable to test stability\n");
          return_value = 0;
          return;
        }                     
        if (its == 10 || its == 20) {
          t += x;
          for (i=1;i<=nn;i++) a[i][i] -= x;
          s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
          y=x=0.75*s;
          w = -0.4375*s*s;
        }
        ++its;
        for (m=(nn-2);m>=l;m--) {
          z=a[m][m];
          r=x-z;
          s=y-z;
          p=(r*s-w)/a[m+1][m]+a[m][m+1];
          q=a[m+1][m+1]-z-r-s;
          r=a[m+2][m+1];
          s=fabs(p)+fabs(q)+fabs(r);
          p /= s;
          q /= s;
          r /= s;
          if (m == l) break;
          u=fabs(a[m][m-1])*(fabs(q)+fabs(r));
          v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
          if ((float)(u+v) == v) break;
        }
        for (i=m+2;i<=nn;i++) {
          a[i][i-2]=0.0;
          if  (i != (m+2)) a[i][i-3]=0.0;
        }
        for (k=m;k<=nn-1;k++) {
          if (k != m) {
            p=a[k][k-1];
            q=a[k+1][k-1];
            r=0.0;
            if (k != (nn-1)) r=a[k+2][k-1];
            if (x=fabs(p)+fabs(q)+fabs(r)) {
              p /= x;
              q /= x;
              r /= x;
            }
          }
          if (s=SIGN(sqrt(p*p+q*q+r*r),p)) {
            if (k == m) {
              if (l != m)
            a[k][k-1] = -a[k][k-1];
            } else
              a[k][k-1] = -s*x;
            p += s;
            x=p/s;
            y=q/s;
            z=r/s;
            q /= p;
            r /= p;
            for (j=k;j<=nn;j++) {
              p=a[k][j]+q*a[k+1][j];
              if (k != (nn-1)) {
            p += r*a[k+2][j];
            a[k+2][j] -= p*z;
              }
              a[k+1][j] -= p*y;
              a[k][j] -= p*x;
            }
            mmin = nn<k+3 ? nn : k+3;
            for (i=l;i<=mmin;i++) {
              p=x*a[i][k]+y*a[i][k+1];
              if (k != (nn-1)) {
            p += z*a[i][k+2];
            a[i][k+2] -= p*r;
              }
              a[i][k+1] -= p*q;
              a[i][k] -= p;
            }
          }
        }
          }
        }
      } while (l < nn-1);
    }
}



           /*------------------------------------*
            *       Active_bits(int i,int p)     *
            *                                    *
            * Returns the number of active bits  *
            * of population 'i' (i.e. the number *
            * of patterns to which pop. 'i' is   *
            * selective).                        *
            * p = total number of patterns.      *
            *------------------------------------*/

int Active_bits(int i, int p)
{
  int k,j,n=0;

  k=1;
  for (j=0;j<p;j++)
    {
      if (i&k) n++;
      k = k<<1;
    }
  return n;
}


           /*------------------------------------*
            *       float gammln(float xx)       *
            *                                    *
            * Returns the value of ln(Gamma(xx)) *
            * for xx > 0                         *
            *          ( By Num. Rec. )          *
            *                                    *
            *------------------------------------*/

#include <math.h>

float gammln(float xx)
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
            -1.231739516,0.120858003e-2,-0.536382e-5};
  int j;
  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}


           /*------------------------------------*
            *         float factln(int n)        *
            *                                    *
            *           Returns ln(n!)           *
            *           Uses gammln();           *
            *          ( By Num. Rec. )          *
            *                                    *
            *------------------------------------*/

#include <stdio.h>

float factln(int n)
{
  static float a[101];
  float gammln();

  if (n < 0) printf("\n>>> Negative factorial in routine factln()\n");
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
  else return gammln(n+1.0);
}


           /*------------------------------------*
            *       float bico(int n,int k)      *
            *                                    *
            * Returns the binomial coefficient   *
            * (n / k) as a floating point number *
            * Uses factln();                     *
            *          ( By Num. Rec. )          *
            *                                    *
            *------------------------------------*/

#include <math.h>

float bico(int n,int k)
{
  float factln();
  return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}
