#pragma once

#include "state.h"
#include <fparticle\particle.h>


/* losses gives the energy loss for particle b(E)=dE/dt */ 
double losses(double E, double r, Particle& particle, State& state);

double adiabaticLosses(double E, double z, Particle& particle);

/* escapeTime gives (t-1_adv + t-1_dec)^-1 */ 
double escapeTime(double E, Particle& particle);

double decayTime(double E, Particle& particle);

double diffusionRate(double E, double r, Particle& particle);



double accelerationRate(double E, Particle& particle);