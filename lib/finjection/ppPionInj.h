#pragma once

#include <fmath/physics.h>
#include <fparticle/Particle.h>



/*pion injection due to pp interaction*/ 
double ppPionInj(double E, const Particle& creator,
	const double density, const SpaceCoord& psc);
	