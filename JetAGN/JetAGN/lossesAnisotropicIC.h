#pragma once

#include <fparticle\particle.h>

#include <string>
#include <iostream>
#include <fstream>


/* lossesAnisotropicIC calculates the energy loss [erg/s] due to anisotropic IC*/
double lossesAnisotropicIC(double E, Particle& particle, double r);
