#pragma once

#include <fparticle/Particle.h>
#include <fmath/physics.h>

#include <iostream>

/*erg/s/cm^3 */

double luminosityIC_old(double E, const Particle& creator, const SpaceCoord& distCoord, fun1 tpf, double phEmin, double phEmax);

double luminosityIC(double E, const Particle& creator, const SpaceCoord& distCoord, const ParamSpaceValues& tpf,
					double tpEmin, double tpEmax);

/*
// Moderski et al. (2005):
double luminosityIC_2(double E, const Particle& creator, const SpaceCoord& distCoord, const ParamSpaceValues& tpf,
					double tpEmin, double tpEmax);
					
double luminositySyKN(double Eph, const Particle& creator, const SpaceCoord& distCoord, double magf);*/