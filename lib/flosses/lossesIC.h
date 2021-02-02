#pragma once


#include <fparticle/Particle.h>




double lossesIC_old(double E, Particle& particle, fun1 tpf, double phEmin, double phEmax);

double lossesIC(double E, Particle& particle, const ParamSpaceValues& tpf, 
				const SpaceCoord& psc, double phEmin, double phEmax);
double lossesIC_Th(double E, Particle& particle, const ParamSpaceValues& tpf, 
				const SpaceCoord& psc, double phEmin, double phEmax);
				
// Moderski et al. (2005):
//double lossesIC(double E, Particle& particle, const ParamSpaceValues& tpf, 
//				const SpaceCoord& psc, double phEmin, double phEmax);