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

	double E = 1.0e9*1.6e-12;  //elijo 1GeV que es la E a la cual las perdidas rad son importantes
	//con este valor reproduzco el 5x10^40 erg/s del paper de valenti
	// electronMass*cLight2;// P2(electronMass*cLight2) / (boltzmann*parameters.starT);

	double Sj = pi*P2(jetRadius(z, openingAngle));
	double stagPoint = stagnationPoint(z);
	double So = 100.0*pi*P2(stagPoint);

	double vel_lat = cLight*openingAngle;

	double tad = adiabaticLosses(E, z, vel_lat) / E; //ad esta en erg/s

//<<<<<<< HEAD
	double wph = Lj *P2(pi) / (100.0*cLight*Sj);//  (cLight*4.0*Sj);
	//double wph1 = parameters.Lj /  (cLight*4.0*Sj);

	double wmag = P2(computeMagField(z)) / (8.0*pi); // [av] ver si se justifica intentar usar el cubo magf, creo que no.
	double tsin = 4.0 * thomson*cLight*wmag*(E / P2(electronMass*cLight2)) / 3.0;
	double tic = 4.0*cLight*thomson*wph*E / P2(electronMass*cLight2);
	double trad = tsin + tic;

	double frad = 1.0 / (1.0 + tad / trad);

	if (So > Sj)
	{
		So = Sj;
	}
	double Lnt = accEfficiency*Lj*(So / Sj)*frad*pow(Dlorentz, 4) / P2(Gamma);
//=======
//	double wph = Lj / (cLight*4.0*Sj);
//	double trad = 4.0*cLight*thomson*wph*E / P2(electronMass*cLight2);
//
//	double frad = 1.0 / (1.0 + tad / trad);
//
//	double Lnt = accEfficiency*Lj*(So / Sj)*frad*pow(Dlorentz, 4) / P2(Gamma);
//>>>>>>> globals

	double f = starDensity(z)*Lnt; //saque el Sj

	return f;   

}

double nonThermalLuminosity(double intRmin, double intRmax)
{
	double integral = intCilindric(intRmin, intRmax, dLnt); // (double z){ return dLnt(z, dummie) });
	return integral;
	//return 1.0e38;
}