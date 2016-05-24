#include "modelParameters.h"

#include "State.h"

#include "lossesAnisotropicIC.h"

#include <fmath\interpolation.h>
#include <fmath\RungeKutta.h>

#include <fparameters\parameters.h>
#include <fmath\physics.h>
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
	return fmagneticField(z, computeModelB0(parameters.Lj, parameters.openingAngle));
}

double jetRadius(double z, double openingAngle)
{
	return z*openingAngle;
}

double eEmax(double z, double B)
{
	double Reff = 10.0*stagnationPoint(z);
	double vel_lat = cLight*parameters.openingAngle;

	double Emax_ad = parameters.accEfficiency*3.0*jetRadius(z, parameters.openingAngle)*cLight*electronCharge*B / (vel_lat*parameters.Gamma);
	double Emax_syn = electronMass*cLight2*sqrt(parameters.accEfficiency*6.0*pi*electronCharge / (thomson*B));
	double Emax_hillas = electronCharge*B*Reff;
	double min1 = std::min(Emax_syn, Emax_syn);

	//std::ofstream file;
	//file.open("Emax.txt", std::ios::out);

	//std::cout << z << '\t' << Emax_ad << '\t' << Emax_syn << '\t' << Emax_hillas << '\t' << std::endl;

	return std::min(min1, Emax_hillas);
		
}

double stagnationPoint(double z)
{
	double Mdot_wind = 1.0e-8*solarMass / yr;
	double v_wind = 2.0e7;

	double stagPoint = sqrt(Mdot_wind*v_wind*cLight / (4.0*parameters.Lj))*jetRadius(z,parameters.openingAngle);

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

void setParameters(boost::property_tree::ptree& cfg)
{
	
	//double mBH = 1.0e7*solarMass;  //black hole mass
	//double rg = mBH*gravitationalConstant / cLight2;

	//double z0 = 50.0*rg; //50 * Rschw  //OJO! si cambian, cambiar tmb nonthermalLuminosity!!

	parameters.starT = cfg.get<double>("starT",3.0e3);

	parameters.Lj = cfg.get<double>("Lj", 1.0e43);
	parameters.openingAngle = cfg.get<double>("openingAngle", 0.1);  //jet opening angle

//	parameters.B0 = sqrt(8.0*parameters.Lj / cLight) / parameters.openingAngle;  //ojo que esto es Bo*z0

//	Rsp = 1.0e14; //distance to stagnation point

	parameters.Gamma = cfg.get<double>("Gamma", 10);

	double inc = 10.0*pi / 180; //ang obs del jet
	double beta = 1.0 - 1.0 / P2(parameters.Gamma);
	parameters.Dlorentz = 1.0 / (parameters.Gamma*(1.0 - cos(inc)*beta));

	parameters.accEfficiency = cfg.get<double>("accEfficiency", 0.1);

	//magneticField = fmagneticField(z0,B0);  //el primero lo calculo en r = z0
	//density = nWindDensity(Rc, starR); //el primero lo calculo en r = Rc

// Data of electrons and protons
	parameters.primaryIndex = cfg.get<double>("primaryIndex", 2.0);
	//factor_qrel   = 3.0; 

//	const auto nEnergyPoints = cfg.get<int>("particle.default.dim.energy.samples", 10);
//
//	ParticleCfg<Electron>::config.logEmin = cfg.get<double>("particle.electron.dim.energy.min", 6.0);
//	ParticleCfg<Electron>::config.logEmax = cfg.get<double>("particle.electron.dim.energy.max", 15.0);
//	ParticleCfg<Electron>::config.nE = cfg.get<int>("particle.electron.dim.energy.samples", nEnergyPoints);
//
////Data of photons
//
//	ParticleCfg<Photon>::config.logEmin = cfg.get<double>("particle.photon.dim.energy.min", -6.0);
//	ParticleCfg<Photon>::config.logEmax = cfg.get<double>("particle.photon.dim.energy.max", 12.0);
//	ParticleCfg<Photon>::config.nE = cfg.get<int>("particle.photon.dim.energy.samples", nEnergyPoints);

	//parameters.targetPhotonEmin = pow(10.0, parameters.photonLogEmin)*1.6e-12;  //0.15e3*1.6e-12;  //photonEmin = 0.15 KeV 
	//parameters.targetPhotonEmax = pow(10.0, parameters.photonLogEmax)*1.6e-12;  //150.0e3*1.6e-12;   //cutEnergy  = 150 KeV

	//parameters.rmin = 1.0*pc;
	//parameters.rmax = 1.0e3*pc;
	//parameters.nR = 10;

	//los parametros de t los comento porque el vector t(i) lo construyo como los crossing times de las celdas xi
//	timeMin = 1.0e-2; 
//	timeMax = (rmax / cLight); // 1.0e11; // rmax / cLight;
//	nTimes = 50;

	DefOpt_IntTriple.samples_x = cfg.get<int>("integrate-3.samples.x", DefOpt_IntTriple.samples_x);
	DefOpt_IntTriple.samples_t = cfg.get<int>("integrate-3.samples.t", DefOpt_IntTriple.samples_t);
	DefOpt_IntTriple.samples_y = cfg.get<int>("integrate-3.samples.y", DefOpt_IntTriple.samples_y);

	DefOpt_RungeKutta.samples_x = cfg.get<int>("runge-kutta-2.samples.x", DefOpt_RungeKutta.samples_x);
	DefOpt_RungeKutta.samples_y = cfg.get<int>("runge-kutta-2.samples.y", DefOpt_RungeKutta.samples_y);

	DefOpt_RungeKuttaSimple.samples_x = cfg.get<int>("runge-kutta-1.samples.x", DefOpt_RungeKuttaSimple.samples_x);
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