#pragma once

#include <fparticle\particle.h>

/* Total loss rate of photons by pair production, Pe'er & Waxman 2004*/ 
double photoPairProduction(double E, const Particle& photon, fun1 tpf, double tpEmin, double tpEmax);