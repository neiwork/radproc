#pragma once

#include <fparticle/Particle.h>

/*neutron injection due to photoHadronic interaction
(Atoyan & Dermer 2003) */

//double neutronPgamma(double E, Vector Nproton, Particle& particle, Particle& proton, fun1 tpf);
double neutronPgamma(double E, Particle& neutron, Particle& proton, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax) ;



