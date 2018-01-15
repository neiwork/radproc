#pragma once

#include <fparticle\particle.h>


/* luminosityAnisotropicIC returns the luminosity erg/s (for [Ne] = erg^-1)
in the jet or blob frame, due to anisotropic IC against a disk.

from Dermer & Schlickneiser ApJ 461, 1993
*/

double luminosityAnisotropicIC(double E, Particle& particle, double z,
	double gamma, fun3 tpf, double starT, const SpaceCoord& psc);


