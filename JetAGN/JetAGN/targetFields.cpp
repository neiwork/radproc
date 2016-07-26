#include "targetFields.h"

#include "modelParameters.h"

#include <fmath\RungeKutta.h>
#include <fmath\physics.h>

#include <fparameters/parameters.h>
#include <fparameters\SpaceIterator.h>
#include <boost/property_tree/ptree.hpp>


double starDensity(double z)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double h_d = GlobalConfig.get<double>("h_d")*pc;
	static const double Nrg = GlobalConfig.get<double>("N_s");
	static const double pseda = GlobalConfig.get<double>("pseda");

	double mBH = 1.0e7*solarMass;  //black hole mass
	double rg = mBH*gravitationalConstant / cLight2;

	double z0 = 50.0*rg;
	double RintMax = h_d;
	double RintMin = z0;

	double vol_frac = P2(openingAngle)/4.0;

	double int_z = pow(RintMax, (3.0 - pseda))/ (3.0 - pseda);

	double Astar = Nrg / (4.0 * pi*int_z); 

	double n_s =  Astar / pow(z, pseda);
	
	return n_s;
}


double starBlackBody(double Ep, double z)  //Cyg A y Mrk
{
	static const double starT = GlobalConfig.get<double>("starT");
	static const double eph_s = GlobalConfig.get<double>("eph_s");
	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double E = Ep*Dlorentz;  //esta seria E_lab

	double tcross = h_d / cLight;
	double wph = eph_s*tcross;

	double lumS = eph_s*pi*P3(h_d)/solarLuminosity;

	double Epeak = boltzmann*starT;
	double Ephmin = Epeak / 1000.0;
	double Ephmax = Epeak*1000.0;

	static const double int_E = RungeKuttaSimple(Ephmin, Ephmax, [&Ephmin,&Epeak](double E){
		return (P3(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
	});  
	
	double normalizacion = wph / int_E; 

	double nph = normalizacion*(P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));


	if (z > h_d){
		return nph*P2(Dlorentz)*P2(h_d / z);
	}
	else{ return nph*P2(Dlorentz); }
}

double starIR(double Ep, double z)  //Cyg A y Mrk
{
	static const double starT = GlobalConfig.get<double>("IRstarT");
	static const double lumIR = GlobalConfig.get<double>("LIR");
	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	static const double R_d = GlobalConfig.get<double>("R_d") *pc;
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double E = Ep*Dlorentz;  //esta seria E_lab

	//double tcross = h_d / cLight;
	
	//IR field

	double wph = lumIR*solarLuminosity / (cLight*pi*P2(R_d)); //lumIR*tcross/vol

	double Epeak = boltzmann*starT;
	double Ephmin = Epeak / 100.0;
	double Ephmax = Epeak*100.0;

	static const double int_E = RungeKuttaSimple(Ephmin, Ephmax, [&Ephmin, &Epeak](double E){
		return (P3(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
	});

	double normalizacion = wph / int_E;

	double nph = normalizacion*(P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));

	if (z > h_d){
		return nph*P2(Dlorentz)*P2(h_d / z);
	}
	else{ return nph*P2(Dlorentz); }

}
	

//double starBlackBody(double Ep, double z)  //M87
//{
//
//	static const double starT = GlobalConfig.get<double>("starT");
//	static const double eph_s = GlobalConfig.get<double>("eph_s");
//	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
//	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
//	
//	double E = Ep*Dlorentz;  //esta seria E_lab
//
//	double wph = (eph_s*solarLuminosity)*starDensity(z)*z / cLight;  //para las gigantes rojas de M87 
//
//	double Epeak = boltzmann*starT;
//	double Ephmin = Epeak / 100.0;
//	double Ephmax = Epeak*100.0;
//
//	static const double int_E = RungeKuttaSimple(Ephmin, Ephmax, [&](double E){
//		return (P3(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
//	});
//
//	double normalizacion = wph / int_E;
//
//	double nph = normalizacion*(P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
//		
//	return nph*P2(Dlorentz);
//}


void tpfPSV(ParamSpaceValues& psv, fun2 tpf, Particle photon, double Lorentz)
{
	//ParamSpaceValues psv(photon.ps);

	photon.ps.iterate([&](const SpaceIterator &i){

		const double E = i.val(DIM_E);
		const double r = i.val(DIM_R);
		
		double Eeval = E*Lorentz;

		double result = tpf(Eeval, r);

		psv.set(i, result);

	} );
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



//void targetPhotonEnergies(double& EphminS, double& EphminCMB)
//{
//	static const double starT = GlobalConfig.get<double>("starT", 3.0e3);
//	static const double Tcmb = GlobalConfig.get<double>("targetPhotonEnergies-Tcmb", 2.7);
//
//	double EpS = boltzmann*starT;
//	EphminS = EpS / 100.0;
//	//double EphmaxS = EpS*Dlorentz * 100.0;
//
//	double EpCMB = boltzmann*Tcmb;
//	EphminCMB = EpCMB / 100.0;
//	//double EphmaxCMB = EpCMB*Dlorentz * 100.0;
//
//}



//double gaussian(double x, double mu, double sigma)
//{
//	//here sigma es sigma^2 of the normal distribution
//
//	double factor = sqrt(2.0*pi*sigma);
//
//	return exp(-P2(x - mu) / (2.0*sigma)) / (factor*x);
//
//}