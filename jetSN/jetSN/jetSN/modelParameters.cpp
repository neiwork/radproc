#include "modelParameters.h"

#include "State.h"

//#include "lossesAnisotropicIC.h"

#include <fmath\interpolation.h>
#include <fmath\RungeKutta.h>

#include <fparameters\parameters.h>
#include <fmath\physics.h>
#include <fmath\configure.h>
#include <iostream>
#include <algorithm>


double stagnationPoint(double z)
{
	static const double E_0 = GlobalConfig.get<double>("E_0");
	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	double Rs = pow((E_0*cLight*pi*P2(openingAngle*z) / (4.0*Lj)), 1.0 / 3.0);

	return Rs;

}


inline double computeModelB0(double Lj, double openingAngle) {
	return sqrt(4.0*Lj / cLight) / openingAngle;  //ojo que esto es Bo*z0
}

inline double fmagneticField(double z, double B_o)
{
	static const double subEq = GlobalConfig.get<double>("subEq");
	return sqrt(subEq)*B_o / z;  //la equiparicion es respecto a Lj ~B^2
}

double computeMagField(double z) {
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double Gamma = GlobalConfig.get<double>("Gamma");

	double Blab = fmagneticField(z, computeModelB0(Lj, openingAngle));
	return Blab/ Gamma;  //este es el B en el sistema del jet
	//return Blab * (5.0e16/z);  //este es el B para m=2

} //la densidad de energia (~B^2) transforma con Gamma^2, 
  //B transforma como Gamma

double jetRadius(double z, double openingAngle)
{
	return z*openingAngle;
}

double eEmax(double z, double B)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");

	double Reff = stagnationPoint(z);
	double vel_lat = cLight*openingAngle;

	double Emax_ad = accEfficiency*3.0*jetRadius(z, openingAngle)*cLight*electronCharge*B / (vel_lat*Gamma); //
	double Emax_syn = electronMass*cLight2*sqrt(accEfficiency*6.0*pi*electronCharge / (thomson*B));
	
	double ampl = Gamma; //factor de amplificaci'on de B en la zona del choque
	double Emax_diff = electronCharge*B*Reff*sqrt(3.0*accEfficiency*ampl/2.0);
	double min1 = std::min(Emax_syn, Emax_ad);
	double min2 = std::min(min1, Emax_diff);

	//std::ofstream file;
	//file.open("Emax.txt", std::ios::out);

	//file << z << '\t' << Emax_ad << '\t' << Emax_syn  << '\t' << std::endl;

	return min2;
		
}



//void derive_parameters_r(double E, double z, double t)
//{
//	double B0{ computeModelB0(parameters.Lj, parameters.openingAngle) };
//	//parameters.radius = jetRadius(z, parameters.openingAngle);
//	parameters.magneticField = fmagneticField(z, B0);
//	//electronLogEmax = log10(eEmax(z, magneticField));
//
//	//Rsp = stagnationPoint(z);
//}

double computeDlorentz(double gamma) {

	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree

	double Dlorentz;
	double beta = sqrt(1.0 - 1.0 / P2(gamma));

	if (inc < openingAngle) {
		Dlorentz = pow((14.0*pow(gamma,4) / 3.0) , (1.0 / 4.0));
	}
	else {
		Dlorentz = 1.0 / (gamma*(1.0 - cos(inc)*beta));
	}
	return Dlorentz;
}


void prepareGlobalCfg()
{
	static const double z_int = GlobalConfig.get<double>("z_int")*pc;
	static const double Gamma = GlobalConfig.get<double>("Gamma");

	GlobalConfig.put("Dlorentz", GlobalConfig.get<double>("Dlorentz", computeDlorentz(Gamma)));
	GlobalConfig.put("Bfield", GlobalConfig.get<double>("Bfield", computeMagField(z_int)));

//	DefOpt_IntLosses.samples_x = GlobalConfig.get<int>("integrate-losses.samples.x", DefOpt_IntLosses.samples_x);
	//DefOpt_IntLosses.samples_t = GlobalConfig.get<int>("integrate-losses.samples.t", DefOpt_IntLosses.samples_t);
	//DefOpt_IntLosses.samples_y = GlobalConfig.get<int>("integrate-losses.samples.y", DefOpt_IntLosses.samples_y);

	fmath_configure(GlobalConfig);
}


//


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

void initializeRPoints(Vector& v, double Rmin, double Rmax)
{

	double R_int = pow((Rmax / Rmin), (1.0 / (v.size() - 1)));

	v[0] = Rmin;

	for (size_t i = 1; i < v.size(); ++i){
		v[i] = v[i - 1] * R_int;
	}

}

void initializeCrossingTimePoints(Vector& time, double rMin, double rMax)
{
	double R_int = pow((rMax / rMin), (1.0 / (time.size()-1)));

	Vector v(time.size()+1, 0.0);

	v[0] = rMin;

	for (size_t i = 0; i < time.size(); ++i){
		v[i+1] = v[i] * R_int;
		double delta_x = v[i+1]-v[i];
		time[i] = delta_x / cLight;  //construyo los t_i como los crossing time de las celdas i
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