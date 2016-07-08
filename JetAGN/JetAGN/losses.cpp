#include "losses.h"

#include "targetFields.h"
//#include <flosses\dataLosses.h>
#include <fparameters\parameters.h>
#include <flosses\lossesSyn.h>
#include "lossesAnisotropicIC.h"
#include <flosses\nonThermalLosses.h>
#include <flosses\lossesIC.h>
//#include <flosses\lossesBrem.h>

#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <map>

double losses(double E, double r, Particle& p, State& st, const SpaceCoord& i)
{	
	static const double starT = GlobalConfig.get<double>("starT");
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	double phEmin = boltzmann*starT*1.0e-2; 
	double phEmax = boltzmann*starT*1.0e2;

	double vel_lat = cLight*openingAngle;

	//double r = i->par.R;
	double B = st.magf.get(i); // parameters.magneticField;
	if(p.id == "electron")	{
		return  lossesSyn(E, B, p)
			+adiabaticLosses(E, r, cLight); 
			+lossesIC(E, p, 
			[&E, &r](double E){
			return starBlackBody(E, r); }, phEmin, phEmax)/P2(Dlorentz);
				//lossesAnisotropicIC(E, p, r);
	}

	throw 0;
}

