#pragma once

#include "state.h"
#include <fparticle\particle.h>


/* losses gives the energy loss for particle b(E)=dE/dt */ 
double losses(double E, Particle& particle, State& state, const SpaceCoord& i);
