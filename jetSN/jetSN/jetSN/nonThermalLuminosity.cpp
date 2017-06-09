#include "nonThermalLuminosity.h"

#include "dynamics.h"
#include "targetFields.h"
#include "modelParameters.h"
//#include "cilindricIntegral.h"

#include "write.h"

#include <flosses\nonThermalLosses.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>
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

	double tad_flow = adiabaticLosses(E, z, vel_lat, Gamma) / E; //VER
	double tad = tad_flow/Gamma; //en el lab

	//<<<<<<< HEAD
	double Rs = stagnationPoint(z);
	double wph2 = 3.0*E_0 / (4.0*P3(Rs));
		
	double tic = 4.0*cLight*thomson* wph2*(E / P2(electronMass*cLight2)) / 3.0;
	
	double trad = tic;

	double frad = 1.0 / (1.0 + tad / trad);
	return frad;

}


double dLnt(double z, double Gc, double z_int)
{
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	static const double Gj = GlobalConfig.get<double>("Gamma");
	//static const double Lj = GlobalConfig.get<double>("Lj");


	double Rs = stagnationPoint(z_int); 

	double y = z / z_int;

	double g = Gc / Gj;

	//double Lsc = 4.0*accEfficiency*cLight*P2(Gj)*jetRamPress(z)*Fe(g,y)*pi*P2(Rs);
	
	double cte = jetRamPress(z_int)*pi*cLight*P2(Rs); // = So/Sj
	double F_e = Fe(g, y);
	double Lsc = accEfficiency*F_e*cte*P2(Gj / P2(Gc))/4.0;

	return Lsc;
}



/*
double nonThermalLuminosity(Particle& p, Vector& Gc)
{
	//static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double theta = GlobalConfig.get<double>("openingAngle");
	
	//static const double starT = GlobalConfig.get<double>("starT");

	//double E = P2(electronMass*cLight2) / (boltzmann*starT) / Gamma; //ver si no a Gc

	const double RMIN = p.ps[DIM_R].first();
	const double RMAX = p.ps[DIM_R].last();
	const int N_R = p.ps[DIM_R].size() - 1;

	double R_int = pow((RMAX / RMIN), (1.0 / N_R));

	double L1 = 0.0;

	double z_int = RMIN;

	std::ofstream file;
	file.open("Lnt.txt", std::ios::out);

	double z = RMIN;

	for (int z_ix = 0; z_ix < N_R; z_ix++) {

		double dz = z*(R_int - 1);

		double beta = sqrt(1.0 - 1.0 / P2(Gc[z_ix]));

		double Dlorentz = computeDlorentz(Gc[z_ix]); // 1.0 / (Gc[z_ix] * (1.0 - cos(inc)*beta));
		double boost = pow(Dlorentz, 4) / P2(Gc[z_ix]);

		//L1 = L1 + dLnt(z,Gc[z_ix],z_int)*frad(E, z)*(pi*P2(jetRadius(z, theta)))*dz;
		double Q = dLnt(z, Gc[z_ix], z_int)*boost;

		L1 = L1 + Q;

		file << z/cLight/yr << '\t' << Q << std::endl;

		z = z*R_int;

	}

	file << "Lnt total" << '\t' << L1 << std::endl;

	std::cout << "Lnt total" << '\t' << L1 << std::endl;

	//double Lnt_total = L1;

	//double integral = dLnt(z_int,Gc)*frad(E, z_int); 
//	double boost = pow(Dlorentz, 4) / P2(Gamma);
//	double Lnt_total = integral*boost;


	//file << "E" << '\t' << E / 1.6e-12 << '\t' << "Lnt total" << '\t' << Lnt_total << std::endl;
	
	file.close();
	
	return L1;
			

}*/