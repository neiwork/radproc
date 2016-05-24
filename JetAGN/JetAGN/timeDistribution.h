#pragma once

#include "state.h"
#include <fparticle\particle.h>

using namespace std; 


double effectiveE(double Ee, double Emax, double t, double r, Particle& p, State& st, const SpaceCoord& i );

double timeDistribution(double Ee, double r, double t, Particle& p, State& st, double Eeff, const SpaceCoord& i );

//void timeDistribution(Particle& p, State& st, int z, double& Eeff);


//void matrixChange(Matrix& M, Matrix& aux, int b_i, int b_j);