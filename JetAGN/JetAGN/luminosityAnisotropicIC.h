#pragma once

#include <fparticle\particle.h>

#include <string>
#include <iostream>
#include <fstream>


/* luminosityAnisotropicIC returns the luminosity erg/cm^3/s due to anisotropic IC*/
double luminosityAnisotropicIC(double E, Particle& particle, double r);


