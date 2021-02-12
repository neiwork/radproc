#include "modelParameters.h"

#include "State.h"


#include <fmath\interpolation.h>
#include <fmath\RungeKutta.h>

#include <fparameters\parameters.h>
#include <fmath\physics.h>
#include <fmath\configure.h>
#include <iostream>
#include <algorithm>

//const DimensionCoord
//DIM_E = 0;

inline double computeModelB0(double Lj, double openingAngle) {
	return sqrt(8.0*Lj / cLight) / openingAngle;  //ojo que esto es Bo*z0
}

inline double fmagneticField(double z, double B_o)
{
	static const double subEq = GlobalConfig.get<double>("subEq");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	double B = B_o / z / (Gamma); 

	return B*sqrt(subEq);
}

double computeMagField(double z) {
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Lj = GlobalConfig.get<double>("Lj");

	return fmagneticField(z, computeModelB0(Lj, openingAngle));
}

double computeDens(double z) {
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Lj = GlobalConfig.get<double>("Lj"); 
	static const double Gamma = GlobalConfig.get<double>("Gamma");

	double n = Lj / (pi*P2(z*openingAngle)* (Gamma - 1.0)*electronMass*cLight*cLight2);
	return n;
}


double jetRadius(double z, double openingAngle)
{
	return z*openingAngle;
}

double eEmax(double mass, double z, double B)
{
	//VER agregar tcross = size/cLight;
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double accEfficiency = GlobalConfig.get<double>("etaEfficiency");

	double size = openingAngle*z;

	double vel_lat = cLight*openingAngle;

	double Emax_ad = accEfficiency*3.0*jetRadius(z, openingAngle)*cLight*electronCharge*B / (vel_lat*Gamma);
	double Emax_syn = mass*cLight2*sqrt(accEfficiency*6.0*pi*electronCharge / (thomson*B));
//	double Emax_hillas = electronCharge*B*size;
	double min1 = std::min(Emax_syn, Emax_ad);

	double Emax_diff = electronCharge*B*size*sqrt(3.0*accEfficiency / 2.0);
	double min2 = std::min(min1, Emax_diff);
	
	double Emin_eV = log10(min2 / 1.6e-12);
	return min2;
			
}



double computeDlorentz(double gamma) {

	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree

																		 //double angle = std::max(openingAngle, inc);
	double Dlorentz;
	double beta = sqrt(1.0 - 1.0 / P2(gamma));

	if (inc < openingAngle) {
		Dlorentz = pow((14.0*pow(gamma, 4) / 3.0), (1.0 / 4.0));
	}
	else {
		Dlorentz = 1.0 / (gamma*(1.0 - cos(inc)*beta));
	}
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

	static const double Rdiss_rg = GlobalConfig.get<double>("Rdiss_rg"); //height of interaction in rg units


	GlobalConfig.put("Dlorentz", computeDlorentz(Gamma));

	GlobalConfig.put("Rdiss", computezInt(Mbh, Rdiss_rg));

	static const double Rdiss = GlobalConfig.get<double>("Rdiss");

	GlobalConfig.put("density", computeDens(Rdiss));

	
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