#include <math.h>
#include "nr.h"
#include "nrutil.h"
#define ALF 1.0e-4               // Ensures sufficient decrease in function value.
#define TOLX 1.0e-20              // Convergence criterion on delta x.
#define float double

void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
	    float *f, float stpmax, int *check, float (*func)(float []))
/* Given an n -dimensional point xold[1..n] , the value of the function and gradient there, fold
   and g[1..n] , and a direction p[1..n] , finds a new point x[1..n] along the direction p from
   xold where the function func has decreased “sufficiently.” The new function value is returned
   in f . stpmax is an input quantity that limits the length of the steps so that you do not try to
   evaluate the function in regions where it is undefined or subject to overflow. p is usually the
   Newton direction. The output quantity check is false (0) on a normal exit. It is true (1) when
   x is too close to xold . In a minimization algorithm, this usually signals convergence and can
   be ignored. However, in a zero-finding algorithm the calling program should check whether the
   convergence is spurious. Some “difficult” problems may require double precision in this routine. */
{
  int i;
  float a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,
    test,tmplam;
  
  *check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=1;i<=n;i++) p[i] *= stpmax/sum;        // Scale if attempted step is too big.
  for (slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];
  if (slope >= 0.0) nrerror("Roundoff problem in lnsrch.");
  test=0.0;                        // compute lambdaMin.
  for (i=1;i<=n;i++) {
    temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;
  alam=1.0;                                 // Always try full Newton step first.
  for (;;) {                                // Start of iteration loop.
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
    *f=(*func)(x);
    if (alam < alamin) {                  // Convergence on deltaX. For zero finding,
      for (i=1;i<=n;i++) x[i]=xold[i];    // the calling program should verify the convergence.
      *check=1;
      return;
    } else if (*f <= fold+ALF*alam*slope) return;         // Sufficient function decrease.
    else {                                                // Backtrack.
      if (alam == 1.0)
	tmplam = -slope/(2.0*(*f-fold-slope));            // First time.
      else {                                              // Subsequent backtracks.
	rhs1 = *f-fold-alam*slope;
	rhs2=f2-fold-alam2*slope;
	a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a == 0.0) tmplam = -slope/(2.0*b);
	else {
	  disc=b*b-3.0*a*slope;
	  if (disc < 0.0) tmplam=0.5*alam;
	  else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
	  else tmplam=-slope/(b+sqrt(disc));
	}
	if (tmplam > 0.5*alam)
	  tmplam=0.5*alam;              // lambda <= 0.5*lambda1
      }
    }
    alam2=alam;
    f2 = *f;
    alam=FMAX(tmplam,0.1*alam);         // lambda >= 0.1*lambda1
  }                                     // Try again.
}
