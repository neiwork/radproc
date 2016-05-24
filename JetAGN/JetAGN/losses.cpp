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

double losses(double E, double r, Particle& p, State& st, const SpaceCoord& i)
{	
	//double r = i->par.R;
	double B = st.magf.get(i); // parameters.magneticField;
	if(p.id == "electron")	{
		return  lossesSyn(E, B, p) +
				adiabaticLosses(E, r, cLight) +
				lossesAnisotropicIC(E, p, r);
	}
//	case PT_proton:
//		return  lossesSyn(E, particle) + lossesHadronics(E, particle)
//			+ lossesPhotoHadronic(E,particle,state.tpf);
//		break;
//	}
	throw 0;
}

