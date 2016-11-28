#include "targetFields.h"

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

	double vol_frac = P2(openingAngle)/4.0;

	double int_z = pow(RintMax, (3.0 - pseda))/ (3.0 - pseda);

	double Astar = Nrg / (4.0 * pi*int_z); // (pi*P2(openingAngle)*int_z);

	double n_s =  Astar / pow(z, pseda);
	
	return n_s;// *vol_frac;
}


double starBBStarburst(double Ep, double z)  //Cyg A y Mrk
{
	static const double starT = GlobalConfig.get<double>("starT");
	static const double eph_s = GlobalConfig.get<double>("eph_s");
	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	static const double R_d = GlobalConfig.get<double>("R_d") *pc;
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double E = Ep*Dlorentz;  //esta seria E_lab

	double tcross = h_d / cLight;

	double wph = eph_s*tcross;

	double Epeak = boltzmann*starT;
	double Ephmin = Epeak / 1000.0;
	double Ephmax = Epeak*1000.0;

	static const double int_E = RungeKuttaSimple(Ephmin, Ephmax, [&Ephmin,&Epeak](double E){
		return (P3(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
	});  
	
	double normalizacion = wph / int_E; 

	double nph = normalizacion*(P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));


	//if (z > R_d/2.0) {
	//	return nph*P2(Dlorentz)*P2(R_d / z/2.0);
	//}
	//else {
		return nph*P2(Dlorentz);
	//}

}
	

double sBBM87(double Ep, double z)  //M87
{

	static const double starT = GlobalConfig.get<double>("starT");
	static const double eph_s = GlobalConfig.get<double>("eph_s");
	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	
	double E = Ep*Dlorentz;  //esta seria E_lab

	
	double Epeak = boltzmann*starT;
	double Ephmin = Epeak / 100.0;
	double Ephmax = Epeak*100.0;

	static const double int_E = RungeKuttaSimple(Ephmin, Ephmax, [&](double E){
		return (P3(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
	});
	
	double nph = (P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
		
	double normalization;
	//para M87 tmb extiendo el jet mas alla del radio del bulge
	if (z > h_d) {
		double wph = (eph_s*solarLuminosity)*starDensity(h_d)*h_d / cLight;  //para las gigantes rojas de M87 
		double norm = wph / int_E;
		normalization = norm;// *P2(h_d / z / 2.0);
	}
	else {
		double wph = (eph_s*solarLuminosity)*starDensity(z)*z / cLight;  //para las gigantes rojas de M87 
		double norm = wph / int_E;
		normalization = norm;
	}

	return normalization*nph*P2(Dlorentz);
}

double starBlackBody(double Ep, double z)
{
	static const std::string id = GlobalConfig.get<std::string>("id");

	if (id == "M87") {
		//std::cout << "M87" << std::endl;
		return sBBM87(Ep, z);
	}
	else {
		return starBBStarburst(Ep, z);
	}
}


double starIR(double Ep, double z)
{   //Cyg A y Mrk
	static const double starT = GlobalConfig.get<double>("IRstarT");
	static const double IRlum = GlobalConfig.get<double>("IRlum");
	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	static const double R_d = GlobalConfig.get<double>("R_d") *pc;
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double E = Ep*Dlorentz;  //esta seria E_lab

	double tcross = h_d / cLight;

	//IR field
	double wph = IRlum*solarLuminosity / (cLight*pi*P2(R_d));

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

double gxM87(double Ep, double z)
{   //M87
	static const double starT = GlobalConfig.get<double>("starT");
	static const double lum = GlobalConfig.get<double>("eph_s");
	static const double h_d = GlobalConfig.get<double>("h_d") *pc; //radio bulge
	static const double R_d = GlobalConfig.get<double>("R_d") *pc; //zmax
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double E = Ep*Dlorentz;  //esta seria E_lab

	double tcross = h_d / cLight;

	//star field
	double wph = 3.0*lum*solarLuminosity / (cLight*pi*P2(R_d));

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

double cmbBlackBody(double Ep, double z)
{
	//nph'(E') = nph(E) Potter, W 2012
	//double E = Dlorentz*Ep;

	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double E = Ep*Dlorentz;  //esta seria E_lab

	double Tcmb = 2.7;
	double Epeak = boltzmann*Tcmb;
	double Ephmin = Epeak/ 100.0; //*Dlorentz 
	double Ephmax = Epeak*100.0; //*Dlorentz

	double nph = 4.0*pi*2.0*P2(E)*exp(-Ephmin / E)*exp(-E/Ephmax) / ((P3(planck*cLight))*
		(exp(E / Epeak) - 1));

	return nph*P2(Dlorentz);
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