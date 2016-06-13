#pragma once


#include <fparticle\Particle.h>
#include "State.h"

#include <fparameters\ParamSpaceValues.h>


/*This function generates the targetPhotonField, necessary to luminosityIC*/
void photonTarget(Particle& p, State& st);

void photonDistribution(Particle& p, State& st);