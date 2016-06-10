#include "nonThermalLuminosity.h"

#include "targetFields.h"
#include "modelParameters.h"
#include "cilindricIntegral.h"

#include <flosses\nonThermalLosses.h>
//#include <fmath\RungeKutta.h>
#include <fparameters\parameters.h>

#include <fmath\physics.h>



double dLnt(double z)  //esta es la función que depende del número de estrellas a tiempo t
{
	
	//double Tstar = 3.0e3; 
	double E = 1.0e9*1.6e-12;  //elijo 1GeV que es la E a la cual las perdidas rad son importantes
	//con este valor reproduzco el 5x10^40 erg/s del paper de valenti
	// electronMass*cLight2;// P2(electronMass*cLight2) / (boltzmann*parameters.starT);

	double Sj = pi*P2(jetRadius(z, parameters.openingAngle));
	double stagPoint = stagnationPoint(z);
	double So = 100.0*pi*P2(stagPoint);

	double vel_lat = cLight*parameters.openingAngle;

	double tad = adiabaticLosses(E, z, vel_lat) / E; //ad esta en erg/s

	double wph = parameters.Lj *P2(pi) / (100.0*cLight*Sj);//  (cLight*4.0*Sj);
	//double wph1 = parameters.Lj /  (cLight*4.0*Sj);

	double wmag = P2(parameters.magneticField) / (8.0*pi);
	double tsin = 4.0 * thomson*cLight*wmag*(E / P2(electronMass*cLight2)) / 3.0;
	double tic = 4.0*cLight*thomson*wph*E / P2(electronMass*cLight2);
	double trad = tsin + tic;

	double frad = 1.0 / (1.0 + tad / trad);

	if (So > Sj)
	{
		So = Sj;
	}
	double Lnt = parameters.accEfficiency*parameters.Lj*(So / Sj)*frad*pow(parameters.Dlorentz, 4) / P2(parameters.Gamma);

	double f = starDensity(z)*Lnt; //saque el Sj

	return f;   

}

double nonThermalLuminosity(double intRmin, double intRmax)
{
	double integral = intCilindric(intRmin, intRmax, dLnt); // (double z){ return dLnt(z, dummie) });
	return integral;
	//return 1.0e38;
}