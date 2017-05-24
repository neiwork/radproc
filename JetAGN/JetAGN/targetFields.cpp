#include "targetFields.h"

#include "modelParameters.h"

#include <iostream>
#include <fmath\RungeKutta.h>
#include <fmath\physics.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>

double f(double h,double z)
{
	static const double R_d = GlobalConfig.get<double>("R_d") *pc;
	double res = 0.5*(h - z)*log(P2(R_d / (z - h)) + 1.0) + R_d*atan((h - z) / R_d);
	return res;
}
double wphStarburst(double z)
{
	static const double eph_s = GlobalConfig.get<double>("eph_s");
	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	//static const double R_d = GlobalConfig.get<double>("R_d") *pc;

	//double tcross = h_d / cLight;

	double cte = eph_s/ (2.0*cLight);

	double wph = cte*(f(h_d, z) - f(-h_d, z));

	return wph;
}

double wph(double z)
{
	static const std::string id = GlobalConfig.get<std::string>("id");

	if (id == "M87") {
		//std::cout << "M87" << std::endl;
		return wphM87(z);
	}
	else {
		return wphStarburst(z);
	}
}


double starBlackBody(double Ep, double z)
{

	static const double starT = GlobalConfig.get<double>("starT");
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double E = Ep*Dlorentz;  //esta seria E_lab

	double Epeak = boltzmann*starT;
	double Ephmin = Epeak / 1000.0;
	double Ephmax = Epeak*1000.0;

	static const double int_E = RungeKuttaSimple(Ephmin, Ephmax, [&Ephmin, &Epeak](double E) {
		return (P3(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
	});

	double normalizacion = wph(z) / int_E;

	double nph = normalizacion*(P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));

	return nph*P2(Dlorentz);
	
}


double wphIR(double z)
{

	static const double IRlum = GlobalConfig.get<double>("IRlum");
	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	static const double R_d = GlobalConfig.get<double>("R_d") *pc;

	//IR field
	double vol = 4.0*pi*P2(R_d)*(2.0*h_d);
	double eph_s = IRlum*solarLuminosity / (vol); //[]=erg/s/cm3

	double cte = eph_s / (2.0*cLight);

	double wph = cte*(f(h_d, z) - f(-h_d, z));

	return wph;
}

double starIR(double Ep, double z)
{   //Cyg A y Mrk
	static const double starT = GlobalConfig.get<double>("IRstarT");
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double E = Ep*Dlorentz;  //esta seria E_lab


	double Epeak = boltzmann*starT;
	double Ephmin = Epeak / 100.0;
	double Ephmax = Epeak*100.0;

	static const double int_E = RungeKuttaSimple(Ephmin, Ephmax, [&Ephmin, &Epeak](double E) {
		return (P3(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
	});

	double normalizacion = wphIR(z) / int_E;

	double nph = normalizacion*(P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));


	//if (z > R_d / 2.0) {
	//	return nph*P2(Dlorentz)*P2(R_d / z / 2.0);
	//}
	//else {
		return nph*P2(Dlorentz);
	//}

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