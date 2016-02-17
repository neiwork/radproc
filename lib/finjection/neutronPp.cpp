#include "neutronPp.h"

#include <flosses\dataLosses.h>
#include <fparameters\parameters.h>
#include <flosses\lossesHadronics.h>
#include <fmath\interpolation.h>


double neutronPp(double E, Vector Nproton, Particle& particle, Particle& proton)  
{
	
	double protonDist = proton.dist(E);// interpol(E, proton.energyPoints, Nproton, Nproton.size() - 1);

	DataLosses data;
	data.E = E;
	data.mass = proton.mass;

	double t_1   = lossesHadronics(E, proton)/E;
	
	double emissivity = 0.5*t_1*protonDist;
	
	return emissivity;
}
