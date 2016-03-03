#pragma once

#include <fparticle\particle.h>


/* this is the Synchrotron luminosity without Synchrotron self absorption; erg/s/cm^3*/ 
double luminositySynchrotron(double E, const Particle& c);

/* this is the Synchrotron luminosity with Synchrotron self absorption*/ 
//double luminositySynchrotron_conSSA(double E, const Particle& creator);


/* this is the Synchrotron luminosity produced by secondary pairs */ 
//double luminositySynchrotronSec(double E, const Particle& creator);