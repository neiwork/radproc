#pragma once

#include "state.h"
#include <fparticle\particle.h>



//void distribution(Particle& particle, State& state);
void distribution(Particle& p, State& st, Vector& Gc, Vector& Rc);

//void distWOLosses(Particle& p, State& st);