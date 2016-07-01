#include "nonThermalLuminosity.h"

#include "targetFields.h"
#include "modelParameters.h"
#include "cilindricIntegral.h"

#include <flosses\nonThermalLosses.h>
#include <fparameters\parameters.h>

#include <fmath\physics.h>

#include <iostream>
#include <boost/property_tree/ptree.hpp>



double frad(double E, double z)
{
	//static const double starT = GlobalConfig.get<double>("starT");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double Lj = GlobalConfig.get<double>("Lj");

	//static const double eph_s = GlobalConfig.get<double>("eph_s"); 
	static const double Ls = GlobalConfig.get<double>("eph_s")*solarLuminosity; //M87
	static const double h_d = GlobalConfig.get<double>("h_d") *pc; 

	double tcross = h_d / cLight;


	double wph = Ls*tcross*starDensity(z);  //para las gigantes rojasd de M87 
	//double wph = eph_s*tcross; // Lj / (cLight*4.0*Sj);


	double Sj = pi*P2(jetRadius(z, openingAngle));

	double vel_lat = cLight*openingAngle;

	double tad_flow = adiabaticLosses(E, z, vel_lat) / E; //ad esta en erg/s

	double tad = tad_flow/Gamma; //en el lab

	//<<<<<<< HEAD

	double wph2 = Lj / (cLight*4.0*Sj);
	
	double wmag = P2(computeMagField(z)*P2(Gamma)) / (8.0*pi); // [av] ver si se justifica intentar usar el cubo magf, creo que no.
	double tsin = 4.0 * thomson*cLight*wmag*(E / P2(electronMass*cLight2)) / 3.0;
	double tic = 4.0*cLight*thomson*wph*(E / P2(electronMass*cLight2)) / 3.0;
	double trad = tsin + tic;

	double frad = 1.0 / (1.0 + tad / trad);
	return frad;

}


double dLnt(double z)  //esta es para M87
{	

	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");

	static const double MdotWind = GlobalConfig.get<double>("Mdot")*solarMass / yr;
	static const double vWind = GlobalConfig.get<double>("vWind");

	//static const double Lj = GlobalConfig.get<double>("Lj");

	//double Sj = pi*P2(jetRadius(z, openingAngle));
	//double stagPoint = stagnationPoint(z);
	//double So = 100.0*pi*P2(stagPoint);

	double ratioS = 100.0*MdotWind*vWind*cLight / (4.0);  ///Lj; == So / Sj si tomo las estrellas iguales
	double Lnt = accEfficiency*(ratioS); //*Lj
	double nstar = starDensity(z);
	double f = nstar*Lnt; //saque el Sj, el doppler boosting y el frad

	return f;

}

/*
double dLnt(double z)  //esta la uso para CygA y Mrk
{

	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	
	static const double int_m = GlobalConfig.get<double>("I");  //erg cm-3 /v

	static const double h_d = GlobalConfig.get<double>("h_d")*pc;  //pc


	double f = accEfficiency*int_m*100.0*cLight*pi / 4.0;

	if (z > h_d){ return 0.0; }
	else{ return f; }

}
*/

double nonThermalLuminosity(double intRmin, double intRmax)
{
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	static const double starT = GlobalConfig.get<double>("starT");

	double E = P2(electronMass*cLight2) / (boltzmann*starT);

	//double E = 1.0e9*1.6e-12;

	double integral = intCilindric(intRmin, intRmax, 
		[&E](double z){return dLnt(z)*frad(E, z); }); // (double z){ return dLnt(z, dummie) });  
	
	double Lnt_total = integral*pow(Dlorentz, 4) / P2(Gamma);

	double Lnt_total_pri = Lnt_total / Gamma;

	std::cout << "Lnt total" << '\t' << Lnt_total << std::endl;
	std::cout << "Lnt FF total" << '\t' << Lnt_total_pri << std::endl;

	return Lnt_total;
		
		

}