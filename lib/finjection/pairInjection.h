#pragma once

#include <fparticle\particle.h>


/* pairInjection calculates the pair injection due to photon-photon annihilation 
Aharonian, Vila & Aha 2009 */ 
double pairInjection(double E, Vector Nphoton, Particle& particle, Particle& photon, fun1 tpf);

/*este es el que uso para las acopladas*/ 
double pairInjectionFlor(double E, Particle& particle, Particle& photon, fun1 tpf);
