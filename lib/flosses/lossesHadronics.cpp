#include "lossesHadronics.h"

#include "crossSectionInel.h"

#include <fparameters/parameters.h>

#include <fmath/physics.h>



double lossesHadronics(double E, double density, Particle& particle)
{
	double inelasticity = 0.5;

	if (particle.id == "pion") {
		return (2.0 / 3.0)*cLight*density*inelasticity*crossSectionHadronic(E)*E;  //mb
	}
	else { //particle == proton or neutron
		return cLight*density*inelasticity*crossSectionHadronic(E)*E;  //mb
	}
}
	
/*
double lossesHadronics(double E, Particle& particle)
{
	double inelasticity = 0.5; 

	switch (particle.type)	{
		case PT_proton:
		case PT_neutron:
			return cLight*density*inelasticity*crossSectionHadronic(E)*E;  //mb
			break;
		case PT_pion:
			return (2.0/3.0)*cLight*density*inelasticity*crossSectionHadronic(E)*E;  //mb
			break;
	}
	throw 0;
}

//el E multiplicando es porque son perdidas y no cooling time

*/