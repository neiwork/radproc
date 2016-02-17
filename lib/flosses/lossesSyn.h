#pragma once

#include <fparticle\particle.h>

class DataLosses;

/* lossesSyn calculates the energy loss due to Synchrotron radiation*/ 
double lossesSyn(double E, Particle& particle );

/* lossesSyn calculates the energy loss due to Synchrotron radiation for secondary pairs*/ 
 double lossesSynSec(double E, Particle& particle);
