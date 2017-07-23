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


double wphStarburst(double z)
{
	static const double eph_s = GlobalConfig.get<double>("eph_s");
	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	static const double R_d = GlobalConfig.get<double>("R_d") *pc;

	double cte = eph_s / (2.0*cLight);

	double wph = cte*(f(h_d, z, R_d) - f(-h_d, z, R_d));

	return wph;
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

	double normalizacion = wphStarburst(z) / int_E;

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

	double normalizacion = wphIR(z) / int_E;

	double nph = normalizacion*(P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));


	//if (z > R_d / 2.0) {
	//	return nph*P2(Dlorentz)*P2(R_d / z / 2.0);
	//}
	//else {
		return nph*P2(Dlorentz);
	//}

}









