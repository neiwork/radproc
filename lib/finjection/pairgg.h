#pragma once

#include <fparticle/Particle.h>


/* pairInjection calculates the pair injection due to photon-photon annihilation 
Aharonian, F. A.; Atoian, A. M.; Nagapetian, A. M. (1983); Vieyro & Romero 2012 */ 

double pairGammaGamma(double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf,
						const SpaceCoord& distCoord, double tpEminSoft, double tpEmaxSoft,
							double tpEminG, double tpEmaxG);
double pairGammaGammaNew(double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf,
						const SpaceCoord& distCoord, double tpEminSoft, double tpEmaxSoft,
							double tpEminG, double tpEmaxG);
double pairGammaGamma2(double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax);