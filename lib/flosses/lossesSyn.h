#pragma once

#include <fparticle/Particle.h>



/* lossesSyn calculates the energy loss [erg/s] due to Synchrotron radiation*/ 
double lossesSyn(double E, double magneticField, Particle& particle);

/* lossesSyn calculates the energy loss due to Synchrotron radiation for secondary pairs*/ 
double lossesSynSec(double E, double magneticField, Particle& particle);
