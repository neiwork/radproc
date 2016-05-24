#include "targetFields.h"


#include "modelParameters.h"
#include <fmath\physics.h>

#include <fparameters/parameters.h>

double starDensity(double z)
{
	static const double openingAngle = GCFG.get<double>("openingAngle", 0.1);

	double Nrg = 4.0e7;

	double pseda = 2.0; // 2.0;

	double mBH = 1.0e7*solarMass;  //black hole mass
	double rg = mBH*gravitationalConstant / cLight2;

	double z0 = 50.0*rg;
	double RintMax = 1.0e3*pc;
	double RintMin = z0;// rmin;// z0;

	double int_z = pow(RintMax, (3.0 - pseda)) - pow(RintMin, (3.0 - pseda)) / (3.0 - pseda);

	double Astar = Nrg / (pi*P2(openingAngle)*int_z);

	double n_s = Astar / pow(z, pseda);
	return n_s;
}


double starBlackBody(double E, double r)
{

	static const double starT = GCFG.get<double>("starT", 3.0e3);
	static const double starR = GCFG.get<double>("starR", 1.0e13);

	double Epeak = boltzmann*starT;
	double Ephmin = Epeak / 100.0;
	//double Ephmax = Epeak*10.0;

	//8.0*pi/(P3(planck*cLight)*P2(starR / r); reemplazo esta constante por otra normalizacion
	//double normalizacion = starDensity(r) / P3(Dlorentz*boltzmann*starT);
	
	double normalizacion = 2.0*P2(starR / r) / P3(planck*cLight);

	double nph =  (P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1))
		*normalizacion;
	
	return nph ; /// P2(Dlorentz)  
}


double cmbBlackBody(double E)
{
	//nph'(E') = nph(E) Potter, W 2012
	//double E = Dlorentz*Ep;

	double Tcmb = 2.7;
	double Epeak = boltzmann*Tcmb;
	double Ephmin = Epeak/ 100.0; //*Dlorentz 
	double Ephmax = Epeak*100.0; //*Dlorentz

	double nph = 4.0*pi*2.0*P2(E)*exp(-Ephmin / E)*exp(-E/Ephmax) / ((P3(planck*cLight))*
		(exp(E / Epeak) - 1));

	return nph ;  
}


void targetPhotonEnergies(double& EphminS, double& EphminCMB)
{
	static const double starT = GCFG.get<double>("starT", 3.0e3);
	static const double Tcmb = GCFG.get<double>("targetPhotonEnergies-Tcmb", 2.7);

	double EpS = boltzmann*starT;
	EphminS = EpS / 100.0;
	//double EphmaxS = EpS*Dlorentz * 100.0;

	double EpCMB = boltzmann*Tcmb;
	EphminCMB = EpCMB / 100.0;
	//double EphmaxCMB = EpCMB*Dlorentz * 100.0;

}