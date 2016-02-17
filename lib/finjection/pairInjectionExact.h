#pragma once

#include <fparticle\particle.h>


/* Gives the electron/positron pair injection 1/erg s cm^3, 
given in Boetcher & Schlickeiser 1997, Peer & Waxman 2004*/ 
double pairInjectionExact(double E, Particle& particle, Particle& photon, fun1 tpf);
