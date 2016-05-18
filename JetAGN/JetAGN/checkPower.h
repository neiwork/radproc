#pragma once


#include "State.h"



/*the following function checks if int( N(E,r)*2*pi*Rj(z)^2 dz ) is below the inyected power*/


double computeInjectedPower(const ParamSpaceValues& dist, int t_ix);
