#pragma once

#include <fparticle/Particle.h>


/* Gives the electron/positron pair injection 1/erg s cm^3, 
given in Boetcher & Schlickeiser 1997, Peer & Waxman 2004*/ 
//double pairInjectionExact(double E, Particle& particle, Particle& photon, fun1 tpf);
double pairInjectionExact(double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax);

//double pairInjection2(double E, Vector Nphoton, Particle& particle, Particle& photon, fun1 tpf);