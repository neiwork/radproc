#pragma once
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

double integrator_qag(gsl_function *, double, double, double, double, size_t, double *,int *);
double integrator_cquad(gsl_function *,double,double,double,double,size_t,double *,int *);
double integrator_qags(gsl_function *,double,double,double,double,size_t,double *,int *);