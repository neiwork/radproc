#include "nonThermalLuminosity.h"

#include "targetFields.h"
#include "modelParameters.h"
#include "cilindricIntegral.h"

#include <flosses\nonThermalLosses.h>
#include <fparameters\parameters.h>

#include <fmath\physics.h>


#include <boost/property_tree/ptree.hpp>


double frad(double E, double z)
{
	//static const double starT = GlobalConfig.get<double>("starT");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double Lj = GlobalConfig.get<double>("Lj");

	double Sj = pi*P2(jetRadius(z, openingAngle));

	double vel_lat = cLight*openingAngle;

	double tad_flow = adiabaticLosses(E, z, vel_lat) / E; //ad esta en erg/s

	double tad = tad_flow*Gamma; //en el lab

	//<<<<<<< HEAD
	double wph = Lj / (cLight*4.0*Sj);

	double wmag = P2(computeMagField(z)*P2(Gamma)) / (8.0*pi); // [av] ver si se justifica intentar usar el cubo magf, creo que no.
	double tsin = 4.0 * thomson*cLight*wmag*(E / P2(electronMass*cLight2)) / 3.0;
	double tic = 4.0*cLight*thomson*wph*(E / P2(electronMass*cLight2)) / 3.0;
	double trad = tsin + tic;

	double frad = 1.0 / (1.0 + tad / trad);
	return frad;

}


double dLnt(double z)  //esta es la función que depende del número de estrellas a tiempo t
{
	

	static const double MdotWind = GlobalConfig.get<double>("Mdot")*solarMass / yr;
	static const double vWind = GlobalConfig.get<double>("vWind");

	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double E = 1.0e6*1.6e-12;  //elijo 1GeV que es la E a la cual las perdidas rad son importantes
	//con este valor reproduzco el 5x10^40 erg/s del paper de valenti
	// electronMass*cLight2;// P2(electronMass*cLight2) / (boltzmann*parameters.starT);

	//double Sj = pi*P2(jetRadius(z, openingAngle));
	//double stagPoint = stagnationPoint(z);
	//double So = 100.0*pi*P2(stagPoint);

	double ratioS = 100.0*MdotWind*vWind*cLight / (4.0*Lj);  //== So / Sj si tomo las estrellas iguales

	double Lnt = accEfficiency*Lj*(ratioS)*frad(E,z)*pow(Dlorentz, 4) / P2(Gamma);

	double f = starDensity(z)*Lnt; //saque el Sj

	return f;   

}

double nonThermalLuminosity(double intRmin, double intRmax)
{
	double integral = intCilindric(intRmin, intRmax, dLnt); // (double z){ return dLnt(z, dummie) });
	return integral;

}