#include <stdio.h>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "integrators.h"

double integrator_qag(gsl_function *F, double a, double b, double epsabs, double epsrel,
					size_t limit, double *error, int *status)
{
	double result = 0.0;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(101);
	gsl_set_error_handler_off();
	
	*status = gsl_integration_qag(F,a,b,epsabs,epsrel,limit,6,w,&result,error);
	
	gsl_integration_workspace_free(w);
	return result;
}

double integrator_qags(gsl_function *F, double a, double b, double epsabs, double epsrel,
					size_t limit, double *error, int *status)
{
	double result = 0.0;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(101);
	gsl_set_error_handler_off();
	
	*status = gsl_integration_qags(F,a,b,epsabs,epsrel,limit,w,&result,error);
	
	gsl_integration_workspace_free(w);
	return result;
}

double integrator_cquad(gsl_function *F, double a, double b, double epsabs, double epsrel,
		size_t limit, double *error, int *status)
{
	double result = 0.0;
	gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(101);
	gsl_set_error_handler_off();
	*status = gsl_integration_cquad(F,a,b,epsabs,epsrel,w,&result,error,&limit);
	
	if (*status) {
		std::cout << "ERROR = " << *status << std::endl;
	}
	
	gsl_integration_cquad_workspace_free(w);
	return result;
}