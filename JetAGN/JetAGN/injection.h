#pragma once

#include <fparticle\Particle.h>
//#include "CoupledEqSys.h"
#include "State.h"

//#include "annihilationInjection.h"
//#include "pairAnnihilation.h"


double primaryInjection(double E, Particle& particle);

//void injection(Particle& particle, Particle& creator, fun1 tpf);

void injection(Particle& p, State& st);


//double electronInjection(double E, const CoupledEqSys* previous, CoupledEqSys* current);

//double positronInjection(double E, const CoupledEqSys* previous, CoupledEqSys* current);



