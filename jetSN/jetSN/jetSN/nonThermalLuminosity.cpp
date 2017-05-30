#include "nonThermalLuminosity.h"

#include "targetFields.h"
#include "modelParameters.h"
//#include "cilindricIntegral.h"

#include "write.h"

#include <flosses\nonThermalLosses.h>
#include <fparameters\parameters.h>

#include <fmath\physics.h>

#include <iostream>
#include <boost/property_tree/ptree.hpp>



double frad(double E, double z) //VER
{
	
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double E_0 = GlobalConfig.get<double>("E_0");  //erg
	
	double vel_lat = cLight*openingAngle;

	double tad_flow = adiabaticLosses(E, z, vel_lat) / E; //ad esta en erg/s
	double tad = tad_flow/Gamma; //en el lab

	//<<<<<<< HEAD
	double Rs = stagnationPoint(z);
	double wph2 = 3.0*E_0 / (4.0*P3(Rs));
		
	double tic = 4.0*cLight*thomson* wph2*(E / P2(electronMass*cLight2)) / 3.0;
	
	double trad = tic;

	double frad = 1.0 / (1.0 + tad / trad);
	return frad;

}


double dLnt(double z)  
{
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	
	
	static const double Lj = GlobalConfig.get<double>("Lj"); 
	static const double theta = GlobalConfig.get<double>("openingAngle");
	
	double Rj = jetRadius(z, theta);
	double Rs = stagnationPoint(z); 
	double sigma = P2(Rs / Rj);


	double f = accEfficiency*Lj*sigma;

	return f;
}


double nonThermalLuminosity()
{
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	static const double z_int = GlobalConfig.get<double>("z_int")*pc;  //pc

	static const double starT = GlobalConfig.get<double>("starT");

	double E = P2(electronMass*cLight2) / (boltzmann*starT) / Gamma;
	
	double integral = dLnt(z_int)*frad(E, z_int); 

	double boost = pow(Dlorentz, 4) / P2(Gamma);
	double Lnt_total = integral*boost;

	std::ofstream file;
	file.open("Lnt.txt", std::ios::out);

	file << "E" << '\t' << E / 1.6e-12 << '\t' << "Lnt total" << '\t' << Lnt_total << std::endl;
	file << "boost" << '\t' << boost << '\t' << "Doppler" << '\t' << Dlorentz  << std::endl;

	file.close();

	std::cout << "Lnt total" << '\t' << Lnt_total << std::endl;

	return Lnt_total;
			

}