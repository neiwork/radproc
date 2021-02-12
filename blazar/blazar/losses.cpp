#include "losses.h"

//#include "targetFields.h"
#include <fparameters\parameters.h>
#include <flosses\lossesSyn.h>
//#include "lossesAnisotropicIC.h"
#include <flosses\nonThermalLosses.h>
#include <flosses\lossesIC.h>
#include <flosses\lossesHadronics.h>
#include <flosses\lossesPhotoHadronic.h>

#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <map>

double losses(double E, Particle& p, State& st, const SpaceCoord& i)
{
	//static const std::string id = GlobalConfig.get<std::string>("id");
	
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double zInt = GlobalConfig.get<double>("Rdiss");
	static const double density = GlobalConfig.get<double>("density");
	

	double vel_lat = cLight*openingAngle;

	//double r = i->par.R;
	double B = st.magf.get(i); // parameters.magneticField;
	
	double loss = 0.0;

	if (p.id == "electron"){
		
		double syn = lossesSyn(E, B, p);

		loss = syn 
			//+ lossesIC(E, p,
			//	[&E, &r](double E) {
			//return starIR(E, r); }, phEmin, phEmax);
			+ adiabaticLosses(E, zInt, vel_lat, Gamma);
			
	}
	else if (p.id == "proton") {

		double pp = lossesHadronics(E, density, p);
		double pgamma = lossesPhotoHadronic(E, p, st.tpf, i, st.photon.emin(), st.photon.emax());

		loss = pp + pgamma
			+ adiabaticLosses(E, zInt, vel_lat, Gamma);

		}
	
	return loss;

}

