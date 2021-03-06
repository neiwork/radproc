#pragma once


#include "State.h"



/*this function checks if sum dr int ( E*Q(E,r)*dE) [erg/s] is below the inyected power*/

double computeInjectedPower(const ParamSpaceValues& dist, int t_ix);


/*this function checks if sum dr int ( E*N(E,r)*dE) [erg] is conserved*/

double computeInjectedEnergy(const ParamSpaceValues& dist);
