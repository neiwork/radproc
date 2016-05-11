#include "nonThermalLosses.h"



#include <fmath\physics.h>



double adiabaticLosses(double E, double z, double vel_lat)  //en [erg/s]
{
	double jetRadius = z*parameters.openingAngle;

	return 2.0*parameters.Gamma*(vel_lat*E / (3.0*jetRadius));
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
