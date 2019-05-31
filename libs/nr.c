//
// nr.c    ---   Oct 2000
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


#include <math.h>
#include "nrutil.h"



// nerf(x) = e^(x^2)*(1+erf(x))

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
#undef a2

#define FUNC(x) ((*func)(x))
double trapzd(double (*func)(double), double a, double b, int n)
{
    double x,tnm,sum,del;
    static double s;
    int it,j;

    if (n == 1) {
        return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
    } else {
        for (it=1,j=1;j<n-1;j++) it <<= 1;
        tnm=it;
        del=(b-a)/tnm;
        x=a+0.5*del;
        for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
        s=0.5*(s+(b-a)*sum/tnm);
        return s;
    }
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz%1&0(9p1,1'9. */

// #define EPS 1.0e-6
#define JMAX 20
double qsimp(double (*func)(double), double a, double b, double eps)
{
    double trapzd(double (*func)(double), double a, double b, int n);
    void nrerror(char error_text[]);
    int j;
    double s,st,ost,os;

   // mod:
   if(a==b) return 0.0;
   
    ost = os = -1.0e30;
    for (j=1;j<=JMAX;j++) {
        st=trapzd(func,a,b,j);
        s=(4.0*st-ost)/3.0;
        if (fabs(s-os) < eps*fabs(os)) return s;
        os=s;
        ost=st;
    }
    // nrerror("Too many steps in routine qsimp");
    return -1.0;
}
// #undef EPS
#undef JMAX
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz%1&0(9p1,1'9. */



#define NRANSI
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50
void sort(unsigned long n, float arr[])
{
    unsigned long i,ir=n,j,k,l=1;
    int jstack=0,*istack;
    float a,temp;

    istack=ivector(1,NSTACK);
    for (;;) {
        if (ir-l < M) {
            for (j=l+1;j<=ir;j++) {
                a=arr[j];
                for (i=j-1;i>=1;i--) {
                    if (arr[i] <= a) break;
                    arr[i+1]=arr[i];
                }
                arr[i+1]=a;
            }
            if (jstack == 0) break;
            ir=istack[jstack--];
            l=istack[jstack--];
        } else {
            k=(l+ir) >> 1;
            SWAP(arr[k],arr[l+1])
            if (arr[l+1] > arr[ir]) {
                SWAP(arr[l+1],arr[ir])
            }
            if (arr[l] > arr[ir]) {
                SWAP(arr[l],arr[ir])
            }
            if (arr[l+1] > arr[l]) {
                SWAP(arr[l+1],arr[l])
            }
            i=l+1;
            j=ir;
            a=arr[l];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i],arr[j]);
            }
            arr[l]=arr[j];
            arr[j]=a;
            jstack += 2;
            if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
            if (ir-i+1 >= j-l) {
                istack[jstack]=ir;
                istack[jstack-1]=i;
                ir=j-1;
            } else {
                istack[jstack]=j-1;
                istack[jstack-1]=l;
                l=i;
            }
        }
    }
    free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz%1&0(9p1,1'9. */


#define EPS1 0.001
#define EPS2 1.0e-8
float probks(float alam)
{
    int j;
    float a2,fac=2.0,sum=0.0,term,termbf=0.0;

    a2 = -2.0*alam*alam;
    for (j=1;j<=100;j++) {
        term=fac*exp(a2*j*j);
        sum += term;
        if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
        fac = -fac;
        termbf=fabs(term);
    }
    return 1.0;
}
#undef EPS1
#undef EPS2
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz%1&0(9p1,1'9. */


// ksone: modified (23 Oct 2000) by S Fusi & G La Camera
// to allow for correlations via the parameter 'float corrn':
// corrn: number of correlated points (~tau_corr/dt) 
#define NRANSI
void ksone(float data[], unsigned long n, float corrn, float (*func)(float), float *d,
    float *prob)
{
    float probks(float alam);
    void sort(unsigned long n, float arr[]);
    unsigned long j;
    float dt,en,ff,fn,fo=0.0;

    sort(n,data);
    en=n;
    *d=0.0;
    for (j=1;j<=n;j++) {
        fn=j/en;
        ff=(*func)(data[j]); // ff: theor cum distr
        dt=FMAX(fabs(fo-ff),fabs(fn-ff)); // fn: exp cum distr
        if (dt > *d) *d=dt;
        fo=fn;
      // printf("%4.4f   %4.4f  %4.4f \n",data[j],fn,ff); // debug
    }

   // (23 Oct 2000: correlations) -------
   en = en/corrn;
   // -----------------------------------
   
    en=sqrt(en);
    *prob=probks((en+0.12+0.11/en)*(*d));
   // printf("%>>> arg of probks: %4.4f   d: %4.4f\n",(en+0.12+0.11/en)*(*d),*d); // debug
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz%1&0(9p1,1'9. */
