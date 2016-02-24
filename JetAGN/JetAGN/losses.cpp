#include "losses.h"

#include <flosses\dataLosses.h>
#include <fparameters\parameters.h>
#include <flosses\lossesSyn.h>
#include "lossesAnisotropicIC.h"
#include <flosses\nonThermalLosses.h>
//#include <flosses\lossesIC.h>
//#include <flosses\lossesBrem.h>

#include <iostream>
#include <map>


double losses(double E, Particle& p, State& st)
{
	double r = p.ps[1][0]; //no quiero el 0, quiero el actual
	
	SpaceIterator* i = p.ps.current;

	//p.ps.dimensions[1]->values[i];// por que no??
	//p.ps[1][i];

	//p.ps[1].par->R;

	
	switch (p.type)	{
	case PT_electron:

		//double z = state.electron.  R;// Cómo hago para saber el z ?
		return  lossesSyn(E, p) +
				adiabaticLosses(E, r, cLight) +
				lossesAnisotropicIC(E, st.electron, r);
			break;
	}

	//	case PT_proton:
//		return  lossesSyn(E, particle) + lossesHadronics(E, particle)
//			+ lossesPhotoHadronic(E,particle,state.tpf);
//		break;
//	}
//	throw 0;
}

