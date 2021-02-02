#include "nonThermalLosses.h"

#include <fparameters/parameters.h>

#include <fmath/physics.h>

#include <boost/property_tree/ptree.hpp>

double adiabaticLosses(double E, double z, double vel_lat, double gamma)  //en [erg/s]
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	double jetRadius = z*openingAngle;

	return gamma*2.0*(vel_lat*E / (3.0*jetRadius));  
	//en el sist lab es sin Gamma
	//termina quedando return 2.0*cLight*E / (3.0*z);
}


double diffusionRate(double E, double radius, double magneticField)   //en [s]^-1
{
	double D = cLight*E / (3.0*electronCharge*magneticField);  //cLight=v; D(E) in the Bohm regime

	return (2.0 * D) / P2(radius);
}

double BohmDiffusionCoeff(double E, double B)
{
	double larmorR = E/(electronCharge*B);
	return 1.0/3.0 * larmorR * cLight;
}


double accelerationRate(double E, double magneticField, double accEfficiency) //en [s]^-1
{
	return accEfficiency*cLight*electronCharge*magneticField / E;
}


double escapeRate(double size, double vel) //en [1/s]
{
	//static const double Gamma = GlobalConfig.get<double>("Gamma");

	return vel /size;  //ver si necesito algun gamma para escribirlo en el FF
}



/*
double diffCoeff_g(double g, Particle& p, double height, double B, double rho)
{
	double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
	double q = GlobalConfig.get<double>("nonThermal.injection.SDA.powerSpectrumIndex");
	double kMin = 1.0/height;
	double vA = B / sqrt(4.0*pi*rho);
	double rL = g*p.mass*cLight2/(electronCharge*B);
	return zeda * (cLight*kMin) *gsl_pow_2(vA/cLight) * pow(rL*kMin,q-2) * g*g;
}

double accelerationRateSDA(double E, Particle& p, double B, double height, double rho)
{
	double gamma = E / (p.mass*cLight2);
	return diffCoeff_g(gamma,p,height,B,rho) / (gamma*gamma);
}
*/