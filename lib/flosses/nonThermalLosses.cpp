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


double accelerationRate(double E, double magneticField, double accEfficiency) //en [s]^-1
{
	return accEfficiency*cLight*electronCharge*magneticField / E;
}


double escapeRate(double size, double vel) //en [1/s]
{
	//static const double Gamma = GlobalConfig.get<double>("Gamma");

	return vel /size;  //ver si necesito algun gamma para escribirlo en el FF
}