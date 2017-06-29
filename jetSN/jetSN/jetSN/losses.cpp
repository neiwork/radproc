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
	

	double vel_lat = cLight*openingAngle;

	//double r = i->par.R;
	static const double B{ st.magf.get(i) };
	
	double loss = lossesSyn(E, B, p);
		//+ adiabaticLosses(E, r, vel_lat, gamma);
	 //las perdidas adiabaticas o de escape las considero en escapeRate




	
	return loss;

}




/*if (E > 9.0e10*1.6e-12 && id != "M87") {
static const double starT = GlobalConfig.get<double>("IRstarT");
double phEmin = boltzmann*starT*1.0e-2;
double phEmax = boltzmann*starT*1.0e2;

loss = loss + lossesIC(E, p,
[&E, &r](double E) {
return starIR(E, r); }, phEmin, phEmax);
}*/
