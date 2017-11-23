#include "operations.h"

#include <cmath>


int sgn(double x)
{
	return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

double safeLog10(double x) {
	return (x>0.0) ? log10(x) : log(1.0);
}