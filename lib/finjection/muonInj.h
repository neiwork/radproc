#pragma once

#include <fparticle/Particle.h>

double muonInj(double E, const Particle& p, const Particle& c, const SpaceCoord& psc);
double muonInjNew(double E, const Particle& p, const Particle& c, const SpaceCoord& psc);

//double muonInj(double E, Particle& particle, Particle& pion_creator);