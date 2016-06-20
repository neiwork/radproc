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

inline double computeModelB0(double Lj, double openingAngle) {
	return sqrt(8.0*Lj / cLight) / openingAngle;  //ojo que esto es Bo*z0
}

inline double fmagneticField(double z, double B_o)
{
	return B_o / z;
}

double computeMagField(double z) {
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double Gamma = GlobalConfig.get<double>("Gamma");

	double Blab = fmagneticField(z, computeModelB0(Lj, openingAngle));
	return Blab / P2(Gamma);  //este es el B en el sistema del jet
}

double jetRadius(double z, double openingAngle)
{
	return z*openingAngle;
}

double eEmax(double z, double B)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");

	double Reff = 10.0*stagnationPoint(z);
	double vel_lat = cLight*openingAngle;

	double Emax_ad = accEfficiency*3.0*jetRadius(z, openingAngle)*cLight*electronCharge*B / (vel_lat); //*Gamma
	double Emax_syn = electronMass*cLight2*sqrt(accEfficiency*6.0*pi*electronCharge / (thomson*B));
	//double Emax_hillas = electronCharge*B*Reff;
	double min1 = std::min(Emax_syn, Emax_syn);

	//std::ofstream file;
	//file.open("Emax.txt", std::ios::out);

	//std::cout << z << '\t' << Emax_ad << '\t' << Emax_syn << '\t' << Emax_hillas << '\t' << std::endl;

	//return std::min(min1, Emax_hillas);
	return min1;
		
}

double stagnationPoint(double z)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Lj = GlobalConfig.get<double>("Lj");

	static const double MdotWind = GlobalConfig.get<double>("Mdot")*solarMass / yr;
	static const double vWind = GlobalConfig.get<double>("vWind");


	double stagPoint = sqrt(MdotWind*vWind*cLight / (4.0*Lj))*jetRadius(z,openingAngle);

	return stagPoint;

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
	double inc = 10.0*pi / 180; //ang obs del jet
	double beta = 1.0 - 1.0 / P2(gamma);
	double Dlorentz = 1.0 / (gamma*(1.0 - cos(inc)*beta));
	return Dlorentz;
}

void prepareGlobalCfg()
{
	static const double Gamma = GlobalConfig.get<double>("Gamma", 10);

	GlobalConfig.put("Dlorentz", GlobalConfig.get<double>("Dlorentz", computeDlorentz(Gamma)));

	//parameters.starT = cfg.get<double>("starT", 3.0e3);

	//parameters.Lj = cfg.get<double>("Lj", 1.0e43);
	//parameters.openingAngle = cfg.get<double>("openingAngle", 0.1); //jet opening angle
	//parameters.Gamma = cfg.get<double>("Gamma", 10);
	//parameters.accEfficiency = cfg.get<double>("accEfficiency", 0.1);
	//parameters.primaryIndex = cfg.get<double>("primaryIndex", 2.0);

	//double mBH = 1.0e7*solarMass;  //black hole mass
	//double rg = mBH*gravitationalConstant / cLight2;
	//double z0 = 50.0*rg; //50 * Rschw  //OJO! si cambian, cambiar tmb nonthermalLuminosity!!
	//	parameters.B0 = sqrt(8.0*parameters.Lj / cLight) / parameters.openingAngle;  //ojo que esto es Bo*z0
	//	Rsp = 1.0e14; //distance to stagnation point

	//magneticField = fmagneticField(z0,B0);  //el primero lo calculo en r = z0
	//density = nWindDensity(Rc, starR); //el primero lo calculo en r = Rc

	//factor_qrel   = 3.0; 

	//parameters.targetPhotonEmin = pow(10.0, parameters.photonLogEmin)*1.6e-12;  //0.15e3*1.6e-12;  //photonEmin = 0.15 KeV 
	//parameters.targetPhotonEmax = pow(10.0, parameters.photonLogEmax)*1.6e-12;  //150.0e3*1.6e-12;   //cutEnergy  = 150 KeV

	//parameters.rmin = 1.0*pc;
	//parameters.rmax = 1.0e3*pc;
	//parameters.nR = 10;

	//los parametros de t los comento porque el vector t(i) lo construyo como los crossing times de las celdas xi
//	timeMin = 1.0e-2; 
//	timeMax = (rmax / cLight); // 1.0e11; // rmax / cLight;
//	nTimes = 50;

	DefOpt_IntLosses.samples_x = GlobalConfig.get<int>("integrate-losses.samples.x", DefOpt_IntLosses.samples_x);
	DefOpt_IntLosses.samples_t = GlobalConfig.get<int>("integrate-losses.samples.t", DefOpt_IntLosses.samples_t);
	DefOpt_IntLosses.samples_y = GlobalConfig.get<int>("integrate-losses.samples.y", DefOpt_IntLosses.samples_y);

	fmath_configure(GlobalConfig);
}


//
//void Particle::initializeLinearEnergyPoints( int n )
//{
//	double Emax  = 1.6e-12*pow(10,logEmax);    
//	double Emin  = 1.6e-12*pow(10,logEmin);
//
//	double E_int = (Emax-Emin)/n;
//
//	energyPoints[0] = Emin;
//
//	for (size_t i=1; i < energyPoints.size() ; ++i){  
//		energyPoints[i] = energyPoints[i-1]+E_int;
//	}
//
//}

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