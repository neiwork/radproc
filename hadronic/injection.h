#pragma once

#include <fparticle\Particle.h>

#include "State.h"




/* Returns Q(E,z) in units of 1/erg/s  (it is multiplied by the volume of each celd)*/

void injection(Particle& p, State& st);

void distribution(Particle& p, State& st);



