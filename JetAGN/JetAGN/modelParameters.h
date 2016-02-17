#pragma once

#include <fmath\mathematics.h>
#include <fparticle\Particle.h>


/* define the inital values of the global parameters*/
void setParameters(void );

void initializeRPoints(Vector& v, double Rmin, double Rmax);

/* change_parameters changes the values of some parameters for the iteration on variable r*/
void derive_parameters_r(double E, double R, double T);

double vWind(double r, double starR);

double starBlackBody(double E, double r);

/* coronaLuminosity transforms the termic photon distribution into luminosity. */ 
double coronaLuminosity(double E);

double coronaFlux(double E);

double thermalPF(double E);

double photonPowerLaw(double E);

double diskFlux(double E);

