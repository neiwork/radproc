#include "modelParameters.h"

#include "State.h"

#include "lossesAnisotropicIC.h"

#include <fmath\interpolation.h>
#include <fmath\RungeKutta.h>

#include <fparameters\parameters.h>
#include <fmath\physics.h>
#include <fmath\configure.h>
#include <iostream>
#include <algorithm>


double stagnationPoint(double z)
{
	static const double E0 = GlobalConfig.get<double>("E0")*solarMass/yr;
	static const double vWind = GlobalConfig.get<double>("vWind");
	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	double cte = 3.0*E0*cLight*P2(openingAngle*z) / (10.0*Lj);

	return pow(cte, 1.0 / 3.0);

}





double jetRadius(double z)
{
	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	double ismP = 1.0e-12;

	double z_recol = sqrt(Lj / (ismP*Gj*Gj*pi*cLight));
	
	double radius;

	if (z < z_recol) {
		radius = z*openingAngle;
	}
	else {
		radius = z_recol*openingAngle;
	}

	return radius;
}

inline double fmagneticField(double z)
{
	static const double subEq = GlobalConfig.get<double>("subEq");
	static const double Lj = GlobalConfig.get<double>("Lj");
		
	return sqrt(subEq)*sqrt(4.0*Lj / cLight) / jetRadius(z);

}

double computeMagField(double z) {
	static const double Gamma = GlobalConfig.get<double>("Gamma");

	double Blab = fmagneticField(z);
	return Blab/ Gamma;  //este es el B en el sistema del jet
	//return Blab * (5.0e16/z);  //este es el B para m=2

} //la densidad de energia (~B^2) transforma con Gamma^2, 
  //B transforma como Gamma



double eEmax(double z, double B)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");

	double Reff = 10.0*stagnationPoint(z);
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



double computeDlorentz(double gamma) {

	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree

	//double angle = std::max(openingAngle, inc);
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
	static const double Gamma = GlobalConfig.get<double>("Gamma", 10);

	GlobalConfig.put("Dlorentz", GlobalConfig.get<double>("Dlorentz", computeDlorentz(Gamma)));


	//magneticField = fmagneticField(z0,B0);  //el primero lo calculo en r = z0
	//density = nWindDensity(Rc, starR); //el primero lo calculo en r = Rc

	//factor_qrel   = 3.0; 

	//parameters.targetPhotonEmin = pow(10.0, parameters.photonLogEmin)*1.6e-12;  //0.15e3*1.6e-12;  //photonEmin = 0.15 KeV 
	//parameters.targetPhotonEmax = pow(10.0, parameters.photonLogEmax)*1.6e-12;  //150.0e3*1.6e-12;   //cutEnergy  = 150 KeV

	//parameters.rmin = 1.0*pc;
	//parameters.rmax = 1.0e3*pc;
	//parameters.nR = 10;



	DefOpt_IntLosses.samples_x = GlobalConfig.get<int>("integrate-losses.samples.x", DefOpt_IntLosses.samples_x);
	DefOpt_IntLosses.samples_t = GlobalConfig.get<int>("integrate-losses.samples.t", DefOpt_IntLosses.samples_t);
	DefOpt_IntLosses.samples_y = GlobalConfig.get<int>("integrate-losses.samples.y", DefOpt_IntLosses.samples_y);

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





//void derive_parameters_r(double E, double z, double t)
//{
//	double B0{ computeModelB0(parameters.Lj, parameters.openingAngle) };
//	//parameters.radius = jetRadius(z, parameters.openingAngle);
//	parameters.magneticField = fmagneticField(z, B0);
//	//electronLogEmax = log10(eEmax(z, magneticField));
//
//	//Rsp = stagnationPoint(z);
//}