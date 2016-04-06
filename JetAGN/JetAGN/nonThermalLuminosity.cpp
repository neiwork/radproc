#include "nonThermalLuminosity.h"

#include "modelParameters.h"
#include <flosses\nonThermalLosses.h>
#include <fmath\RungeKutta.h>
#include <fparameters\parameters.h>

#include <fmath\physics.h>


double starDensity(double z)
{
	double Nrg = 4.0e7;
	double n_s = Nrg / P3(z); //star density
	return n_s;
}

double dLnt(double z)  //esta es la función que depende del número de estrellas a tiempo t
{
	
	double Tstar = 3.0e3; //VER
	double E = P2(electronMass*cLight2) / (boltzmann*Tstar);

	double Sj = pi*P2(jetRadius(z, openingAngle));
	double So = 100.0*pi*P2(Rsp);

	double tad = adiabaticLosses(E, z, cLight) / E; //ad esta en erg/s

	double wph = Lj / (cLight*4.0*Sj);
	double trad = 4.0*cLight*thomson*wph*E / P2(electronMass*cLight2);

	double frad = 1.0 / (1.0 + tad / trad);

	double Lnt = accEfficiency*Lj*(So / Sj)*frad*pow(Dlorentz, 4) / P2(Gamma);

	double f = Sj*starDensity(z)*Lnt;

	return f;   //VER ingrese esto de prueba

}

double nonThermalLuminosity ()
{
	double integral = RungeKuttaSimple(rmin, rmax, dLnt); // (double z){ return dLnt(z, dummie) });
	return integral;
	//return 1.0e38;
}