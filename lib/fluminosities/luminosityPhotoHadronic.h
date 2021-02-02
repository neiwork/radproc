#pragma once

#include <fparticle/Particle.h>

/*(Atoyan & Dermer 2003; Vila & Aharonian 2009) */
double luminosityPhotoHadronic(double E, Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& psc,
								double tpEmin, double tpEmax);
