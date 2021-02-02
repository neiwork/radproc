#include <math.h>
#include "nr.h"
#define FACTOR 1.6
#define NTRY 50
//#define float double

//                double (*func)(double)
// std::function< double        (double) > func

int zbrac(std::function<double(double)> func,double *x1,double *x2)
/* Given a function func and an initial guessed range x1 to x2, the routine expands the range
   geometrically until a root is bracketed by the returned values x1 and x2 (in which case zbrac
   returns 1) or until the range becomes unacceptably large (in which case zbrac returns 0). */
{
  void nrerror(char error_text[]);
  int j;
  double f1,f2;
  
  if (*x1 == *x2) nrerror("Bad initial range in zbrac");
  f1=func(*x1);
  f2=func(*x2);
  for (j=1;j<=NTRY;j++) {
    if (f1*f2 < 0.0) return 1;
    if (fabs(f1) < fabs(f2))
      f1=func(*x1 += FACTOR*(*x1-*x2));
    else
      f2=func(*x2 += FACTOR*(*x2-*x1));
  }
  return 0;
}
