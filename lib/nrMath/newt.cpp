#include <math.h>
#include "nr.h"
#include "nrutil.h"
#define MAXITS 3000
#define TOLF 1.0e-20
#define TOLMIN 1.0e-20
#define TOLX 1.0e-20
#define STPMX 10000.0
//#define float double
/* Here MAXITS is the maximum number of iterations; TOLF sets the convergence criterion on
   function values; TOLMIN sets the criterion for deciding whether spurious convergence to a
   minimum of fmin has occurred; TOLX is the convergence criterion on deltax; STPMX is the scaled
   maximum step length allowed in line searches. */


#define FREERETURN {free_dvector(fvec2,1,n);free_dvector(xold,1,n);	\
    free_dvector(p,1,n);free_dvector(g,1,n);free_dmatrix(fjac,1,n,1,n);	\
    free_ivector(indx,1,n);return;}
	
	
int nn2;                   //era nn pero se repetia en brodyn  // Global variables to communicate with fmin.  /nn, fvec nrfunc->agregue un 2 en el nombre
double *fvec2;
//void (*nrfuncv2)(int n, double v[], double f[]);
//void (*nrfuncv)(int n, double v[], double f[]);
std::function<void(int, double [], double [])> nrfuncv2;



void fdjac(int n, double x[], double fvec[], double **df,
	std::function<void(int, double [], double [])> vecfunc);
	
void newt(double x[], int n, int *check,
	  void (*vecfunc)(int, double [], double []))
/* Given an initial guess x[1..n] for a root in n dimensions, find the root by a globally convergent
   Newton’s method. The vector of functions to be zeroed, called fvec2[1..n] in the routine
   below, is returned by the user-supplied routine vecfunc(n,x,fvec2). The output quantity
   check is false (0) on a normal return and true (1) if the routine has converged to a local
   minimum of the function fmin defined below. In this case try restarting from a different initial
   guess. */
{
  //void fdjac(int n, double x[], double fvec2[], double **df,
	//     void (*vecfunc)(int, double [], double []));
  double fmin1(double x[]);
  void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
	      double *f, double stpmax, int *check, double (*func)(double []));
  void lubksb(double **a, int n, int *indx, double b[]);
  void ludcmp(double **a, int n, int *indx, double *d);
  int i,its,j,*indx;
  double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;

  indx=ivector(1,n);
  fjac=dmatrix(1,n,1,n);
  g=dvector(1,n);
  p=dvector(1,n);
  xold=dvector(1,n);
  fvec2=dvector(1,n);                // Define global variables.
  nn2=n;
  nrfuncv2=vecfunc;
  f=fmin1(x);                       // fvec2 is also computed by this call.
  test=0.0;                        // Test for initial guess being a root. Use
  for (i=1;i<=n;i++)               // more stringent test than simply TOLF.
    if (fabs(fvec2[i]) > test) test=fabs(fvec2[i]);
  if (test < 0.01*TOLF) {
    *check=0;
    FREERETURN
  }
  for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);           // Calculate stpmax for line searches.
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
  for (its=1;its<=MAXITS;its++) {                        // Start of iteration loop.
    fdjac(n,x,fvec2,fjac,vecfunc);
    /* If analytic Jacobian is available, you can replace the routine fdjac below with your
       own routine. */
    for (i=1;i<=n;i++) {           // Compute ∇f for the line search.
      for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec2[j];
      g[i]=sum;
    }
    for (i=1;i<=n;i++) xold[i]=x[i];         // Store x,
    fold=f;                                  // and f .
    for (i=1;i<=n;i++) p[i] = -fvec2[i];      // Right-hand side for linear equations.
    ludcmp(fjac,n,indx,&d);                  // Solve linear equations by LU decomposition.
    lubksb(fjac,n,indx,p);
    lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin1);
    // lnsrch returns new x and f.
    //It also calculates fvec2 at the new x when it calls fmin.
    test=0.0;                              // Test for convergence on function values.
    for (i=1;i<=n;i++)
      if (fabs(fvec2[i]) > test) test=fabs(fvec2[i]);
    if (test < TOLF) {
      *check=0;
      FREERETURN
	}
    if (*check) {                   // Check for gradient of f zero, i.e., spurious convergence.
      test=0.0;
      den=FMAX(f,0.5*n);
      for (i=1;i<=n;i++) {
	temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
	if (temp > test) test=temp;
      }
      *check=(test < TOLMIN ? 1 : 0);
      FREERETURN
    }
    test=0.0;                  // Test for convergence on δx.
    for (i=1;i<=n;i++) {
      temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX) FREERETURN
		       }
  nrerror("MAXITS exceeded in newt");
}


#include <math.h>
#include "nrutil.h"
#define EPS 1.0e-4              // Approximate square root of the machine precision.



void fdjac(int n, double x[], double fvec[], double **df,
	std::function<void(int, double [], double [])> vecfunc)
	
//void fdjac(int n, double x[], double fvec2[], double **df,
//	   void (*vecfunc)(int, double [], double []))
/* Computes forward-difference approximation to Jacobian. On input, x[1..n] is the point at
   which the Jacobian is to be evaluated, fvec2[1..n] is the vector of function values at the
   point, and vecfunc(n,x,f) is a user-supplied routine that returns the vector of functions at
   x . On output, df[1..n][1..n] is the Jacobian array. */
{
  int i,j;
  double h,temp,*f;
  f=dvector(1,n);
  for (j=1;j<=n;j++) {
    temp=x[j];
    h=EPS*fabs(temp);
    if (h == 0.0) h=EPS;
    x[j]=temp+h;              // Trick to reduce finite precision error.
    h=x[j]-temp;
	vecfunc(n,x,f);
    //(*vecfunc)(n,x,f);
    x[j]=temp;
    for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec2[i])/h;           // Forward difference formula.
  }
  free_dvector(f,1,n);
}


#include "nrutil.h"


//estas estan definidas arriba?
//extern int nn2;
//extern double *fvec2;
//extern void (*nrfuncv2)(int n, double v[], double f[]);

double fmin1(double x[])
/* Returns f at x. The global pointer *nrfuncv points to a routine that returns the
   vector of functions at x. It is set to point to a user-supplied routine in the calling program.
   Global variables also communicate the function values back to the calling program. */
{
  int i;
  double sum;
 // (*nrfuncv2)(nn2,x,fvec2);
  nrfuncv2(nn2,x,fvec2);
  for (sum=0.0,i=1;i<=nn2;i++) sum += SQR(fvec2[i]);
  return 0.5*sum;
}
