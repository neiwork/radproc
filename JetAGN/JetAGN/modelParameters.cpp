#include "modelParameters.h"

#include "State.h"

#include <fmath\interpolation.h>
#include <fparameters\parameters.h>
#include <fmath\physics.h>
//double thermalPF(double E);


double fmagneticField(double z, double B_o)
{	
	return B_o/z;
}

double jetRadius(double z, double openingAngle)
{
	return z*openingAngle;
}

double eEmax(double z, double B)
{
	double vel = cLight;
	double Emax_ad = accEfficiency*3.0*z*cLight*electronCharge*B / vel;
	double Emax_syn = accEfficiency*6.0*pi*electronCharge*electronMass*cLight2/(thomson*B);
	return std::min(Emax_syn, Emax_syn);
}



void derive_parameters_r(double E, double z, double t)
{
	radius = jetRadius(z, openingAngle);
	magneticField = fmagneticField(z, B0); 
	//electronLogEmax = log10(eEmax(z, magneticField));

	//	density = nWindDensity(r, starR); //el primero lo calculo en r = Rc ?
}


void setParameters(void )
{
	
	double mBH = 14.8*1.99e33;  //black hole mass
	double rg = mBH*gravitationalConstant / cLight2;

	double z0 = 100.0*rg; //50 * Rschw

	double Lj = 1.0e43;
	openingAngle = 0.1;  //jet opening angle

	B0 = sqrt(8.0*Lj / cLight) / openingAngle;  //ojo que esto es Bo*z0

	double inc = 10.0*pi / 180; //ang obs del jet
	Gamma = 10;
	Dlorentz = 1.0 / (Gamma*(1.0 - cos(inc)));

	accEfficiency = 0.1; 

	//magneticField = fmagneticField(z0,B0);  //el primero lo calculo en r = z0
	//density = nWindDensity(Rc, starR); //el primero lo calculo en r = Rc

// Data of electrons and protons
	primaryIndex  = 2.0;
	//factor_qrel   = 3.0; 

	electronLogEmin = 6.0;
	electronLogEmax = 14.0;  


//Data of photons

	photonLogEmin = -2.0;
	photonLogEmax = 12.0;

	targetPhotonEmin = pow(10.0,photonLogEmin)*1.6e-12;  //0.15e3*1.6e-12;  //photonEmin = 0.15 KeV 
	targetPhotonEmax = pow(10.0,photonLogEmax)*1.6e-12;  //150.0e3*1.6e-12;   //cutEnergy  = 150 KeV

	rmin = 1.0*pc;
	rmax = 1.0e3*pc;
	nR = 2;

	timeMin = 1.0e-2;
	timeMax = 1.0e5; // rmax / cLight;
	nTimes = 2;

	nEnergies = 20;        //massive particles
	nPhotonEnergies = 20;  //259;  //photons
}

void initializeRPoints(Vector& v, double Rmin, double Rmax)
{
	//double Emax = 1.6e-12*pow(10, logEmax);
	//double Emin = 1.6e-12*pow(10, logEmin);

	double R_int = pow((Rmax / Rmin), (1.0 / (v.size() - 1)));

	v[0] = Rmin;

	for (size_t i = 1; i < v.size(); ++i){
		v[i] = v[i - 1] * R_int;
	}

}

void initializeLinearPoints(Vector& v, double tMin, double tMax)
{
	double tStep = (tMax - tMin) / (v.size() - 1); // nTimePoints;

	v[0] = tMin;

	for (size_t i=1; i < v.size(); ++i){  
	
		v[i] = v[i-1]+tStep;
	}
}


double blackBody(double E, double r)
{
	double starT = 1.0e5; //VER
	double starR = 1.0e-3*pc;

	double Epeak = boltzmann*starT;
	double Ephmin = Epeak / 100.0;
	//double Ephmax = Epeak*10.0;

	return 8.0*pi*P2(E)*exp(-Ephmin / E) / ((P3(planck*cLight))*
		(exp(E / Epeak) - 1))*P2(starR / r);
}

//double photonPowerLaw(double E)
//{
//	double photonIndex = 2.2;
//	double Aph         = 7.9e8;  // esta es la cte que normaliza la distribucion de fotones termicos
//	double photonEmin  = 0.2e3*1.6e-12;  //photonEmin = 0.15 KeV 
//	double cutEnergy   = 200.0e3*1.6e-12;   //cutEnergy  = 150 KeV 
//
//	return Aph*pow(E,(-photonIndex))*exp(-E/cutEnergy)*exp(-photonEmin/E);
//}
//
//double diskBlackBody(double E)
//{
//	return 8*pi*P2(E)/(P3(planck*cLight)*(exp(E/diskT)-1));	
//}
//
//
//double thermalPF(double E)   //en [erg^-1 cm^-3]
//{
//	return photonPowerLaw(E)+diskBlackBody(E);
//}
//
//
//double starBlackBody (double E, double r)
//{
//	//double E = i.par.E;
//	double Epeak = boltzmann*starT;
//	double Ephmin = Epeak / 100.0;
//	//double Ephmax = Epeak*10.0;
//
//	return 8.0*pi*P2(E)*exp(-Ephmin / E) / ((P3(planck*cLight))*
//		(exp(E / Epeak) - 1))*P2(starR / r);
//}




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