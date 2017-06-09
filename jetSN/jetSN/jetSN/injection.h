#pragma once

#include <fparticle\Particle.h>
//#include "CoupledEqSys.h"
#include "State.h"

//#include "annihilationInjection.h"
//#include "pairAnnihilation.h"



/* Returns Q(E,z) in units of 1/erg/s  (it is multiplied by the volume of each celd)*/

void injection(Particle& p, State& st, Vector& Gc);





