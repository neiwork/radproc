#include "targetFields.h"

#include "modelParameters.h"

#include <iostream>
#include <fmath\RungeKutta.h>
#include <fmath\physics.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>


double f(double h, double z, double R)
{
	double res = 0.5*(h - z)*log(P2(R / (z - h)) + 1.0) + R*atan((h - z) / R);
	if (h == z) { res = 0.0; }
	return res;
}

double eph(const char* id)
{
	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	static const double R_d = GlobalConfig.get<double>("R_d") *pc;

	double starT;
	double eph_s;

	if (id == "IR") {
		starT = GlobalConfig.get<double>("IRstarT");
		static const double IRlum = GlobalConfig.get<double>("IRlum");
		double vol = pi*P2(R_d)*(2.0*h_d);
		eph_s = IRlum*solarLuminosity / (vol); //[]=erg/s/cm3
	}
	else
	{
		starT = GlobalConfig.get<double>("starT");
		eph_s = GlobalConfig.get<double>("eph_s");
	}

	return eph_s;
}


double wphIR(double z, const char* id)
{

	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	static const double R_d = GlobalConfig.get<double>("R_d") *pc;

	double eph_s = eph(id); //[]=erg/s/cm3

	double cte = eph_s / (2.0*cLight);

	double wph = cte*(f(h_d, z, R_d) - f(-h_d, z, R_d));

	return wph;
}


double starBlackBody(double Ep, double z, double gamma)
{

	static const double starT = GlobalConfig.get<double>("starT");
	//static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double Dlorentz = computeDlorentz(gamma);

	double E = Ep*Dlorentz;  //esta seria E_lab

	double Epeak = boltzmann*starT;
	double Ephmin = Epeak / 1000.0;
	double Ephmax = Epeak*1000.0;

	static const double int_E = RungeKuttaSimple(Ephmin, Ephmax, [&Ephmin, &Epeak](double E) {
		return (P3(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
	});

	double normalizacion = wphIR(z, "star") / int_E;

	double nph = normalizacion*(P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));

	return nph*P2(Dlorentz);

}

double starIR(double Ep, double z, double gamma)
{   
	static const double starT = GlobalConfig.get<double>("IRstarT");
	//static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double Dlorentz = computeDlorentz(gamma);

	double E = Ep*Dlorentz;  //esta seria E_lab
	
	double Epeak = boltzmann*starT;
	double Ephmin = Epeak / 100.0;
	double Ephmax = Epeak*100.0;

	static const double int_E = RungeKuttaSimple(Ephmin, Ephmax, [&Ephmin, &Epeak](double E) {
		return (P3(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
	});

	double normalizacion = wphIR(z,"IR") / int_E;

	double nph = normalizacion*(P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));


	//if (z > R_d / 2.0) {
	//	return nph*P2(Dlorentz)*P2(R_d / z / 2.0);
	//}
	//else {
		return nph*P2(Dlorentz);
	//}

}

double nph_ICani2(double E, double gamma, double eta, const char* id)
{
	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	double starT;

	if (id == "IR") {
		starT = GlobalConfig.get<double>("IRstarT");
	}
	else {
		starT = GlobalConfig.get<double>("starT");
	}
	
	double theta = boltzmann*starT / (electronMass*cLight2);
	double theta_min = theta / 100.0;

	double cte = eph(id) / (2.0*cLight);
	
	double K_1 = (electronMass*cLight2)*pow(pi*theta, 4) / 15.0;

	double wph = cte / K_1;
	//esto va por dR, pero ese esta afuera [wph = cm-3]


	double corr = gamma*(1.0 - beta(gamma)*eta);

	double arg = E*corr / theta;
	double den = exp(arg) - 1.0;
	
	if (den > 0.0) {
		return wph*(P2(E) / den);
	}
	else{ return 0.0; }


}

double nph_ICani(double Ep, double z, double r, double gamma, const char* id)
{
	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	double starT;

	if (id == "IR") {
		starT = GlobalConfig.get<double>("IRstarT");
	} else {
		starT = GlobalConfig.get<double>("starT");
	}
	//static const double int_E = RungeKuttaSimple(Ephmin, Ephmax, [&Ephmin, &Epeak](double E) {
	//	return (P3(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
	//});
		
	double cte = eph(id) / (2.0*cLight);
	double theta = boltzmann*starT / (electronMass*cLight2);
	
	double K_1 = (electronMass*cLight2)*pow(pi*theta, 4) / 15.0;

	double wph = cte/K_1;
	//esto va por dR, pero ese esta afuera

	return wph;


}










