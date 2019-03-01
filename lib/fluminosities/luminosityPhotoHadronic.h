#pragma once

#include <fparticle\particle.h>

/*(Atoyan & Dermer 2003; Vila & Aharonian 2009) */
double luminosityPhotoHadronic(double E, const Particle& creator, fun1 tpf, 
	                           const SpaceCoord& psc, double tpEmin, double tpEmax);
