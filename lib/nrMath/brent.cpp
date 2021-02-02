#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "brent.h"

double brent(gsl_function *F, double x_lo, double x_hi, int *status1, int *status2)
{
	gsl_set_error_handler_off();
	double r=0.0;
	int iter = 0, max_iter = 1000;
	const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
	gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

	gsl_root_fsolver_set(s, F, x_lo, x_hi);
	do {
		iter++;
		*status1 = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		*status2 = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);

    } while (*status2 == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);
	return r;
}