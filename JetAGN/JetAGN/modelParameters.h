#pragma once

#include <fmath\mathematics.h>
#include <fparticle\Particle.h>


/* define the inital values of the global parameters*/
void setParameters(void );

void initializeRPoints(Vector& v, double Rmin, double Rmax);

void initializeLinearPoints(Vector& v, double tMin, double tMax);

/* change_parameters changes the values of some parameters for the iteration on variable r*/
void derive_parameters_r(double E, double R, double T);

double blackBody(double E, double r);
double jetRadius(double z, double openingAngle);

//double vWind(double r, double starR);

/* coronaLuminosity transforms the termic photon distribution into luminosity. */ 
double coronaLuminosity(double E);

double coronaFlux(double E);

double thermalPF(double E);

double photonPowerLaw(double E);

double diskFlux(double E);

