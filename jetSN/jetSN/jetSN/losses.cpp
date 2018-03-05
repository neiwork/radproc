#include "losses.h"

#include "targetFields.h"
//#include <flosses\dataLosses.h>
#include <fparameters\parameters.h>
#include <flosses\lossesSyn.h>
//#include "lossesAnisotropicIC.h"
#include <flosses\nonThermalLosses.h>
#include <flosses\lossesIC.h>
//#include <flosses\lossesBrem.h>

#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <map>



double losses(double E, double r, Particle& p, State& st, const SpaceCoord& i, double gamma)
{
	
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gj = GlobalConfig.get<double>("Gamma");

	double vel_lat = cLight*openingAngle;

	double B = computeMagField(r, gamma);

	double loss = lossesSyn(E, B, p)
		+ adiabaticLosses(E, r, vel_lat, gamma); 
	
	return loss;

}

