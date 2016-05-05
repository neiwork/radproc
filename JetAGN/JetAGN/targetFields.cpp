#include "targetFields.h"


#include "modelParameters.h"
#include <fmath\physics.h>



double starDensity(double z)
{
	double Nrg = 4.0e7;

	double pseda = 1.0; // 2.0;

	double RintMax = rmax;

	double Astar = Nrg*(3.0 - pseda) / (pow(RintMax, (3.0 - pseda)) * pi*P2(openingAngle)); //(pow(rmax, (3.0 - pseda)) - pow(rmin, (3.0 - pseda))) / pi;

	double n_s = Astar / pow(z, pseda);
	return n_s;
}

void targetPhotonEnergies(double& EphminS, double& EphminCMB)
{
	double starT = 1.0e4;
	double EpS = boltzmann*starT;
	EphminS = EpS*Dlorentz / 100.0;
	//double EphmaxS = EpS*Dlorentz * 100.0;

	double Tcmb = 2.7;
	double EpCMB = boltzmann*Tcmb;
	EphminCMB = EpCMB*Dlorentz / 100.0;
	//double EphmaxCMB = EpCMB*Dlorentz * 100.0;

}



double starBlackBody(double Ep, double r)
{

	double E = Dlorentz*Ep;

	double starT = 1.0e4;
	double starR = 1.0e-3*pc;

	double Epeak = boltzmann*starT;
	double Ephmin = Epeak / 100.0;
	//double Ephmax = Epeak*10.0;

	//8.0*pi/(P3(planck*cLight)*P2(starR / r); reemplazo esta constante por otra normalizacion
	//double normalizacion = starDensity(r) / P3(Dlorentz*boltzmann*starT);
	
	double normalizacion = 2.0*P2(pc / r) / P3(planck*cLight);

	double nph =  (P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1))
		*normalizacion;
	
	return nph / P2(Dlorentz);  //devuelvo el nph'(E')
}



double cmbBlackBody(double Ep)
{
	//dado que los calculos son en el sistema del jet, llamo a la funcion con Ep

	double E = Dlorentz*Ep;

	double Tcmb = 2.7;
	double Epeak = boltzmann*Tcmb;
	double Ephmin = Epeak*Dlorentz / 100.0;
	//double Ephmax = Epeak*10.0;

	double nph = 4.0*pi*2.0*P2(E)*exp(-Ephmin / E) / ((P3(planck*cLight))*
		(exp(E / Epeak) - 1));

	return nph / P2(Dlorentz);  //devuelvo el nph'(E')
}
