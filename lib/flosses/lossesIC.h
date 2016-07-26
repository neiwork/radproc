#pragma once


#include <fparticle\particle.h>

#include <string>
#include <iostream>
#include <fstream>


//double lossesIC(double E, Particle& particle, fun1 tpf, double phEmin, double phEmax);
double lossesIC(double E, Particle& particle, const SpaceCoord& distCoord, ParamSpaceValues tpf, double phEmin, double phEmax);