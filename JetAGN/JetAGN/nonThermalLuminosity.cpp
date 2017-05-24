#include "nonThermalLuminosity.h"

#include "targetFields.h"
#include "modelParameters.h"
#include "cilindricIntegral.h"

#include "write.h"

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
	//static const double Lj = GlobalConfig.get<double>("Lj");

	//static const double eph_s = GlobalConfig.get<double>("eph_s");
	//static const double h_d = GlobalConfig.get<double>("h_d") *pc;

	static const std::string id = GlobalConfig.get<std::string>("id");

	
	double Sj = pi*P2(jetRadius(z, openingAngle));

	double vel_lat = cLight*openingAngle;

	double tad_flow = adiabaticLosses(E, z, vel_lat) / E; //ad esta en erg/s

	double tad = tad_flow/Gamma; //en el lab

	//<<<<<<< HEAD

	//double wph2 = Lj / (cLight*4.0*Sj);
		
	double tic = 4.0*cLight*thomson* wph(z)*(E / P2(electronMass*cLight2)) / 3.0;
	
	double trad = tic;

	double frad = 1.0 / (1.0 + tad / trad);
	return frad;

}


double dLntM87(double z)  
{	//esta es para M87

	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	static const double h_d = GlobalConfig.get<double>("h_d")*pc;  //pc
	static const double MdotWind = GlobalConfig.get<double>("Mdot")*solarMass / yr;
	static const double vWind = GlobalConfig.get<double>("vWind");

	static const double Lj = GlobalConfig.get<double>("Lj");

	//double Sj = pi*P2(jetRadius(z, openingAngle));
	//double stagPoint = stagnationPoint(z);
	//double So = 100.0*pi*P2(stagPoint);

	double ratioS = 100.0*MdotWind*vWind*cLight / (4.0) / Lj;  // == So / Sj si tomo las estrellas iguales
	double Lnt = accEfficiency*(ratioS)*Lj; //
	double nstar = starDensity(z);
	double f = nstar*Lnt;

	double prueba = ratioS*nstar;

	//return f;
	if (z > h_d)
	{
		return 0.0;
	}
	//else{ return MdotWind*nstar; }
	else { return f; }
}


double dLntStarburst(double z)  
{//esta la uso para CygA y Mrk

	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	
	static const double int_m = GlobalConfig.get<double>("I");  //erg cm-3 /cm
	static const double h_d = GlobalConfig.get<double>("h_d")*pc;  //pc
	static const double Lj = GlobalConfig.get<double>("Lj");  
	
	//double prueba = int_m*100.0*cLight / 4.0/Lj;

	double f = accEfficiency*int_m*100.0*cLight / 4.0;

	if (z > h_d)
	{ return 0.0; }
	//else{ return prueba; }
	else{ return f; }
}

double dLnt(double z)
{
	static const std::string id = GlobalConfig.get<std::string>("id");

	if (id == "M87") {
		//std::cout << "M87" << std::endl;
		return dLntM87(z);
	}
	else {
		return dLntStarburst(z);
	}
}


double nonThermalLuminosity(double intRmin, double intRmax)
{
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	static const double starT = GlobalConfig.get<double>("starT");

	double E = P2(electronMass*cLight2) / (boltzmann*starT) / Gamma;
	
	double integral = intCilindric(intRmin, intRmax,
		[&E](double z){return dLnt(z)*frad(E, z); });
		//[&E](double z){return dLnt(z); });
	
	double boost = pow(Dlorentz, 4) / P2(Gamma);
	double Lnt_total = integral*boost;

	std::ofstream file;
	file.open("Lnt.txt", std::ios::out);

	file << "E" << '\t' << E / 1.6e-12 << '\t' << "Lnt total" << '\t' << Lnt_total << std::endl;
	file << "boost" << '\t' << boost << '\t' << "Doppler" << '\t' << Dlorentz  << std::endl;

	file.close();

	double Lnt_total_pri = Lnt_total / Gamma;

	std::cout << "Lnt total" << '\t' << Lnt_total << std::endl;
	std::cout << "Lnt FF total" << '\t' << Lnt_total_pri << std::endl;

	return Lnt_total;
		
		

}