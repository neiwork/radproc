#pragma once

#include <fparticle\Particle.h>



/* lossesAnisotropicIC calculates the energy loss [erg/s] due to anisotropic IC
due to anisotropic IC against a disk

from Carson & Chiang (2007, astro-ph)
*/


double lossesAnisotropicIC(double E, Particle& particle, double z,
	double gamma, fun3 tpf, double starT);

