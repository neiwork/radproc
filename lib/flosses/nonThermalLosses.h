#pragma once


#include <fparticle\particle.h>

/*vel is the lateral jet expansion*/
double adiabaticLosses(double E, double z, double vel);  //en [erg/s]

/*Diffusion rate in the Bohm regime*/
double diffusionRate(double E, double radius, double magneticField);   //en [s]^-1


double accelerationRate(double E, double magneticField, double accEfficiency); //en [s]^-1