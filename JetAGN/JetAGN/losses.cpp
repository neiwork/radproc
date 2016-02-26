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
	SpaceIterator* i = p.ps.current;
	
	double r = i->par.R;
	double B = i->par.magneticField;
	
	switch (p.type)	{
	case PT_electron:

		return  lossesSyn(E, B, p) +
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

