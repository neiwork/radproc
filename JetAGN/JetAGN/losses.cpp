#include "losses.h"

#include <flosses\dataLosses.h>
#include <fparameters\parameters.h>
#include <flosses\lossesSyn.h>
#include <flosses\nonThermalLosses.h>
//#include <flosses\lossesIC.h>
//#include <flosses\lossesBrem.h>

#include <iostream>
#include <map>


double losses(double E, Particle& particle, State& state)
{
//	switch (particle.type)	{
//	case PT_electron:
	
	//particle.ps.
	
	//double z = state.electron.  R;// Cómo hago para saber el z ?
	//return  lossesSyn(E, particle) + adiabaticLosses(E, z, vel);// +lossesIC(E, particle, state.tpf);
	return 0.0;
//		break;
//	case PT_proton:
//		return  lossesSyn(E, particle) + lossesHadronics(E, particle)
//			+ lossesPhotoHadronic(E,particle,state.tpf);
//		break;
//	}
//	throw 0;
}

