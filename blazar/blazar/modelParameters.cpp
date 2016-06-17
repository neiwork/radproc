#include "modelParameters.h"

#include "State.h"


#include <fmath\interpolation.h>
#include <fmath\RungeKutta.h>

#include <fparameters\parameters.h>
#include <fmath\physics.h>
#include <fmath\configure.h>
#include <iostream>
#include <algorithm>


inline double computeModelB0(double Lj, double openingAngle) {
	return sqrt(8.0*Lj / cLight) / openingAngle;  //ojo que esto es Bo*z0
}

inline double fmagneticField(double z, double B_o)
{
	return B_o / z;
}

double computeMagField(double z) {
	static const double openingAngle = GlobalConfig.get<double>("openingAngle", 0.1);
	static const double Lj = GlobalConfig.get<double>("Lj", 1.0e43);

	return fmagneticField(z, computeModelB0(Lj, openingAngle));
}

double jetRadius(double z, double openingAngle)
{
	return z*openingAngle;
}

double eEmax(double z, double B)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle", 0.1);
	static const double Gamma = GlobalConfig.get<double>("Gamma", 10);
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency", 0.1);

	double size = 0.0; //VER calcular el tama—o
	double vel_lat = cLight*openingAngle;

	double Emax_ad = accEfficiency*3.0*jetRadius(z, openingAngle)*cLight*electronCharge*B / (vel_lat*Gamma);
	double Emax_syn = electronMass*cLight2*sqrt(accEfficiency*6.0*pi*electronCharge / (thomson*B));
	double Emax_hillas = electronCharge*B*size;
	double min1 = std::min(Emax_syn, Emax_syn);


	return std::min(min1, Emax_hillas);
		
}

//double stagnationPoint(double z)
//{
//	static const double openingAngle = GlobalConfig.get<double>("openingAngle", 0.1);
//	static const double Lj = GlobalConfig.get<double>("Lj", 1.0e43);
//
//	double Mdot_wind = 1.0e-8*solarMass / yr;
//	double v_wind = 2.0e7;
//
//	double stagPoint = sqrt(Mdot_wind*v_wind*cLight / (4.0*Lj))*jetRadius(z,openingAngle);
//
//	return stagPoint;
//
//}


double computeDlorentz(double gamma) {
	double inc = 10.0*pi / 180; //ang obs del jet
	double beta = 1.0 - 1.0 / P2(gamma);
	double Dlorentz = 1.0 / (gamma*(1.0 - cos(inc)*beta));
	return Dlorentz;
}

double computezInt(double Mbh, double height) {

	double rg = gravitationalConstant*Mbh*solarMass/cLight2; 
	
	double zInt = height*rg;

	return zInt;
}

void prepareGlobalCfg()
{
	static const double Gamma = GlobalConfig.get<double>("Gamma");

	static const double Mbh = GlobalConfig.get<double>("Mbh");  //black hole mass in solar masses

	static const double height = GlobalConfig.get<double>("height"); //height of interaction in rg units


	GlobalConfig.put("Dlorentz", computeDlorentz(Gamma));

	GlobalConfig.put("zInt", computezInt(Mbh, height));

	
	//DefOpt_IntLosses.samples_x = GlobalConfig.get<int>("integrate-losses.samples.x", DefOpt_IntLosses.samples_x);
	//DefOpt_IntLosses.samples_t = GlobalConfig.get<int>("integrate-losses.samples.t", DefOpt_IntLosses.samples_t);
	//DefOpt_IntLosses.samples_y = GlobalConfig.get<int>("integrate-losses.samples.y", DefOpt_IntLosses.samples_y);

	fmath_configure(GlobalConfig);
}




void initializeEnergyPoints(Vector& v, double logEmin, double logEmax)
{
	double Emax = 1.6e-12*pow(10, logEmax);
	double Emin = 1.6e-12*pow(10, logEmin);

	double E_int = pow((10 * Emax / Emin), (1.0 / (v.size() - 1)));

	v[0] = Emin;

	for (size_t i = 1; i < v.size(); ++i){
		v[i] = v[i - 1] * E_int;
	}

}



/*

double vWind(double r, double starR) //x es la energia
{
	double vInfty = 2.0e8;
	return vInfty*(1.0 - starR / r);
}

double nWindDensity(double r, double starR)
{
	double Mdot = 3.0e-6*solarMass;
	return Mdot / (4.0*pi*P2(r)*vWind(r, starR)*protonMass);
}*/