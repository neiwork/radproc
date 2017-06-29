#include "nonThermalLuminosity.h"

#include "dynamics.h"
//#include "modelParameters.h"

//#include "targetFields.h"
//#include "cilindricIntegral.h"

//#include "write.h"

//#include <flosses\nonThermalLosses.h>
//#include <fparameters\SpaceIterator.h>
//#include <fparameters\Dimension.h>
#include <fparameters\parameters.h>


#include <fmath\physics.h>

#include <iostream>
#include <boost/property_tree/ptree.hpp>


/*
double frad(double E, double z) //VER
{
	
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double E_0 = GlobalConfig.get<double>("E_0");  //erg
	
	double vel_lat = cLight*openingAngle;

	double tad_flow = adiabaticLosses(E, z, vel_lat, Gamma) / E; //VER
	double tad = tad_flow/Gamma; //en el lab

	//<<<<<<< HEAD
	double Rs = blobRadius(z);
	double wph2 = 3.0*E_0 / (4.0*P3(Rs));
		
	double tic = 4.0*cLight*thomson* wph2*(E / P2(electronMass*cLight2)) / 3.0;
	
	double trad = tic;

	double frad = 1.0 / (1.0 + tad / trad);
	return frad;

}*/


double dLnt(double z, double Gc, double z_int, double Rs)
{
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double theta = GlobalConfig.get<double>("openingAngle");


	//double Rs = blobRadius(z_int,z,Gc); 

	double y = z / z_int;

	double g = Gc / Gj;

	//double Lsc = 4.0*accEfficiency*cLight*P2(Gj)*jetRamPress(z)*Fe(g,y)*pi*P2(Rs);
	
	//double cte1 = jetRamPress(z_int)*pi*cLight*P2(Rs); // = So/Sj
	//double F_e = Fe(g, y);
	//double Lsc1 = accEfficiency*F_e*cte1*P2(Gj / P2(Gc))/4.0;
		

	double frac = (Rs/jetRadius(z, theta));

	if (frac > 1.0) {
		Rs = jetRadius(z, theta); 
		std::cout << "SN radius larger than jet radius" << '\t' << 	frac << std::endl;
	}


	double cte = jetRamPress(z)*pi*cLight*P2(Rs); // = So/Sj

	double Lsc = accEfficiency*(1.0/P2(g)-P2(g))*cte/(4.0*P2(Gj));

	return Lsc;
}


