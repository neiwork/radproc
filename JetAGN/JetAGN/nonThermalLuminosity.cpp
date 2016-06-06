#include "nonThermalLuminosity.h"

#include "targetFields.h"
#include "modelParameters.h"
#include "cilindricIntegral.h"

#include <flosses\nonThermalLosses.h>
#include <fparameters\parameters.h>

#include <fmath\physics.h>


#include <boost/property_tree/ptree.hpp>

double dLnt(double z)  //esta es la función que depende del número de estrellas a tiempo t
{
	
	static const double starT = GlobalConfig.get<double>("starT", 3.0e3);
	static const double openingAngle = GlobalConfig.get<double>("openingAngle", 0.1);
	static const double Gamma = GlobalConfig.get<double>("Gamma", 10);
	static const double Lj = GlobalConfig.get<double>("Lj", 1.0e43);
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency", 0.1);
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double E = P2(electronMass*cLight2) / (boltzmann*starT);

	double Sj = pi*P2(jetRadius(z, openingAngle));
	double stagPoint = stagnationPoint(z);
	double So = 100.0*pi*P2(stagPoint);

	double vel_lat = cLight*openingAngle;

	double tad = adiabaticLosses(E, z, vel_lat) / E; //ad esta en erg/s

	double wph = Lj / (cLight*4.0*Sj);
	double trad = 4.0*cLight*thomson*wph*E / P2(electronMass*cLight2);

	double frad = 1.0 / (1.0 + tad / trad);

	double Lnt = accEfficiency*Lj*(So / Sj)*frad*pow(Dlorentz, 4) / P2(Gamma);

	double f = starDensity(z)*Lnt; //saque el Sj

	return f;   

}

double nonThermalLuminosity(double intRmin, double intRmax)
{
	double integral = intCilindric(intRmin, intRmax, dLnt); // (double z){ return dLnt(z, dummie) });
	return integral;
	//return 1.0e38;
}