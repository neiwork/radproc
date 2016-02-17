#pragma once

#include "state.h"
#include <fparticle\particle.h>

using namespace std; 

void timeDistribution(Particle& p, State& st, int z, double& Eeff);
//void timeDistribution(Particle& particle, State& state);

void matrixChange(Matrix& M, Matrix& aux, int b_i, int b_j);