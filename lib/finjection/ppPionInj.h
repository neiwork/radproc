#pragma once

#include <fmath\physics.h>
#include <fparticle\particle.h>



/*pion injection due to pp interaction*/ 
double ppPionInj(double E, const Particle& creator,
	const double density, const SpaceCoord& psc);