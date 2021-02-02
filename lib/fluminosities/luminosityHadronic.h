#pragma once

#include <fparticle\particle.h>


/* Gamma ray emissivity for pp colissions, for a power-law htadron distribution
  From Kelner, Aharonian & Bugayov, 2006.*/

/*[erg s^-1 cm^-3 ]*/
double luminosityHadronic(double E, const Particle& creator,
	const double density, const SpaceCoord& psc);

//double luminosityHadronic(double E, const Particle& creator,
//	const ParamSpaceValues& denf, const SpaceCoord& psc);