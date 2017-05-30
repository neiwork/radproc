#pragma once

#include "state.h"
#include <fparticle\particle.h>


/* losses gives the energy loss for particle b(E)=dE/dt */ 
double losses(double E, double r, Particle& particle, State& state, const SpaceCoord& i);


/* escapeTime gives (t-1_adv + t-1_dec)^-1 */ 
//double escapeTime(double E, Particle& particle);
