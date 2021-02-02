#pragma once

#include <fparticle\particle.h>


/* this is the Synchrotron luminosity without Synchrotron self absorption; [L] = erg^2/s * [N]*/ 
double luminositySynchrotron2(double E, const Particle& c, const SpaceCoord& psc, double magf);

double luminositySynchrotron(double E, const Particle& c, const SpaceCoord& distCoord, const ParamSpaceValues& magf);

/* this is the Synchrotron luminosity with Synchrotron self absorption*/ 
//double luminositySynchrotron_conSSA(double E, const Particle& creator);


/* this is the Synchrotron luminosity produced by secondary pairs */ 
//double luminositySynchrotronSec(double E, const Particle& creator);