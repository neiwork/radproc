#include "targetFieldM87.h"

#include "modelParameters.h"

#include <iostream>
#include <fmath\RungeKutta.h>
#include <fmath\physics.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>




double starDensity(double z) //M87
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double h_d = GlobalConfig.get<double>("h_d")*pc;
	static const double Nrg = GlobalConfig.get<double>("N_s");
	static const double pseda = GlobalConfig.get<double>("pseda");

	double RintMax = h_d;

	//double int_z = pow(RintMax, (3.0 - pseda)) - pow(RintMin, (3.0 - pseda)) / (3.0 - pseda);

	double vol_frac = P2(openingAngle) / 4.0;

	double int_z = pow(RintMax, (3.0 - pseda)) / (3.0 - pseda);

	double Astar = Nrg / (4.0 * pi*int_z); // (pi*P2(openingAngle)*int_z);

	double n_s = Astar / pow(z, pseda);

	return n_s;// *vol_frac;
}


double wph_int(double Rb, double z, double pseda)
{
	//double t_min = 0.0;
	//double t_max = pi;
	//int nt = 100;
	//double dt = (t_max - t_min) / nt;
	double r_min = 1.0e15;  //Rsch 10 ^ 8Msun = 10 ^ 13
	double r_max = Rb;
	int nr = 500;
	double r_int = pow((r_max / r_min), (1.0 / nr));

	double Q = 0.0;

	//double t = t_min;

	//for (size_t i = 0; i < nt; ++i) {

		double r = r_min;

		for (size_t j = 0; j < nr; ++j) {

			double dr = r*(r_int - 1.0);

			//double f = z*pow(r, (2 - pseda))*cos(t)*sin(t)*sqrt(1.0 + P2(r / (z*cos(t))));
			double f = pow(r, (1 - pseda))*log(abs((r + z) / (r - z)));

			Q = Q + f*dr;

			r = r*r_int;

		}

	//	t = t + dt;
	//}

	return Q;
}

double wphM87(double z)
{
	static const double eph_s = GlobalConfig.get<double>("eph_s");
	static const double Rb = GlobalConfig.get<double>("h_d") *pc;
	static const double pseda = GlobalConfig.get<double>("pseda");
	
	double Astar = starDensity(z)*pow(z, pseda);
	double cte = (eph_s*solarLuminosity)*Astar / (2.0*cLight*z);

	double Q = wph_int(Rb, z, pseda);

	return Q*cte;

}


double cmbBlackBody(double Ep, double z)
{
	//nph'(E') = nph(E) Potter, W 2012
	//double E = Dlorentz*Ep;

	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double E = Ep*Dlorentz;  //esta seria E_lab

	double Tcmb = 2.7;
	double Epeak = boltzmann*Tcmb;
	double Ephmin = Epeak / 100.0; //*Dlorentz 
	double Ephmax = Epeak*100.0; //*Dlorentz

	double nph = 4.0*pi*2.0*P2(E)*exp(-Ephmin / E)*exp(-E / Ephmax) / ((P3(planck*cLight))*
		(exp(E / Epeak) - 1));

	return nph*P2(Dlorentz);
}





double gxM87(double Ep, double z)
{   //M87
	static const double starT = GlobalConfig.get<double>("starT");
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	//static const double vWind = GlobalConfig.get<double>("vWind");
	static const double Lj = GlobalConfig.get<double>("Lj");

	double vWind = 2.0e7;
	double Mdot = 1.0e-7*solarMass / yr;

	double E = Ep*Dlorentz;  //esta seria E_lab

	double Rs = (z*0.1)*sqrt(Mdot*vWind*cLight / (4.0*Lj));

	//double tcross = Reff / cLight;
	//double vol = 4.0*pi*P3(Reff) / 3.0;

	double wph = 1.0e2*solarLuminosity / (4.0*cLight*pi*P2(3.0*Rs));

	double Epeak = boltzmann*starT;
	double Ephmin = Epeak / 100.0;
	double Ephmax = Epeak*100.0;

	static const double int_E = RungeKuttaSimple(Ephmin, Ephmax, [&Ephmin, &Epeak](double E) {
		return (P3(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
	});

	double normalizacion = wph / int_E;

	double nph = normalizacion*(P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));

	//if (z > R_d / 2.0) {
	//	return nph*P2(Dlorentz)*P2(R_d / z / 2.0);
	//}
	//else {
	return nph*P2(Dlorentz);
	//}

}





