#pragma once
#include <gsl/gsl_math.h>

double brent(gsl_function *F, double x_lo, double x_hi, int *status1, int *status2);