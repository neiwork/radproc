#pragma once

#include "state.h"
#include <fparticle\particle.h>


/* losses gives the energy loss for particle b(E)=dE/dt */ 
double losses(double E, double r, 
	Particle& p, State& st, const SpaceCoord& i, double gamma);


/* escapeTime gives (t-1_adv + t-1_dec)^-1 */ 
//double escapeTime(double E, Particle& particle);
