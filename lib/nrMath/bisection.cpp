#include "bisection.h"

double bisection (double a, double b, double allerr, int maxmitr, double (*fun)(double))
{
  int itr = 0;
  double x, x1;
  
  while( itr<maxmitr)
    {
      itr++;
      x = (a+b)/2;
      if ((*fun)(a) * (*fun)(x) < 0)
	b=x;
      else
	a=x;
      x1 = (a+b)/2;
      if (fabs(x1-x) < allerr) {
	return x1;
      }
    }
  return x;
}
