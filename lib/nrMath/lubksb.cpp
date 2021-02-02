#include "nr.h"
#define float double

void lubksb(float **a, int n, int *indx, float b[])
/* Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix
   A but rather as its LU decomposition, determined by the routine ludcmp . indx[1..n] is input
   as the permutation vector returned by ludcmp . b[1..n] is input as the right-hand side vector
   B, and returns with the solution vector X. a , n , and indx are not modified by this routine
   and can be left in place for successive calls with different right-hand sides b . This routine takes
   into account the possibility that b will begin with many zero elements, so it is efficient for use
   in matrix inversion. */
{
  int i,ii=0,ip,j;
  float sum;
  
  for (i=1;i<=n;i++) {              // When ii is set to a positive value, it will become the
    ip=indx[i];                     // index of the first nonvanishing element of b. We now
    sum=b[ip];                      // do the forward substitution, equation. The
    b[ip]=b[i];                     // only new wrinkle is to unscramble the permutation as we go.
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;             // A nonzero element was encountered, so from now on we
    b[i]=sum;                       // will have to do the sums in the loop above.
  }
  for (i=n;i>=1;i--) {              // Now we do the backsubstitution.
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];               // Store a component of the solution vector X.
  }                       // All done!
}
