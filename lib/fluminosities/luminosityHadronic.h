#pragma once

#include <fparticle\particle.h>

/*[erg s^-1 cm^-3 ]*/
double luminosityHadronic(double E, const Particle& creator,
	const ParamSpaceValues& denf, const SpaceCoord& psc);