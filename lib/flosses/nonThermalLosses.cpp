#include "nonThermalLosses.h"

#include <fparameters/parameters.h>

#include <fmath/physics.h>

#include <boost/property_tree/ptree.hpp>

double adiabaticLosses(double E, double z, double vel_lat)  //en [erg/s]
{
	static const double openingAngle = GCFG.get<double>("openingAngle", 0.1);
	static const double Gamma = GCFG.get<double>("Gamma", 10);

	double jetRadius = z*openingAngle;

	return 2.0*Gamma*(vel_lat*E / (3.0*jetRadius));
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
