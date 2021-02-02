#include "neutronPp.h"

#include <fparameters\parameters.h>
#include <flosses\lossesHadronics.h>




double neutronPp(double E, Particle& proton ,
	const double density, const SpaceCoord& psc)

{
	
	double protonDist = proton.distribution.interpolate({ { 0, E } }, &psc);  
		
	double t_1   = lossesHadronics(E, density, proton)/E;
	
	double emissivity = 0.5*t_1*protonDist;
	
	return emissivity;
}
