#include "lossesSyn.h"


//#include <gsl/gsl_sf_dilog.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>


double fKN_sy(double b)
{
	double g = (0.5*b + 6.0 + 6.0/b) * log(1.0+b);
	g = g - (11.0/12.0 * b*b*b + 6.0*b*b + 9.0*b + 4.0) / P2(1.0+b);
//	g = g - 2.0 + 2.0 * gsl_sf_dilog(-b);
	return (b > 1e-3 ? 9.0*g/(b*b*b) : 1.0);
}


double lossesSyn(double E, double magf, Particle& p)
{
	double mass = p.mass;
	double gamma = E / (p.mass*cLight2);

	double wmag = magf*magf / (8.0 * pi);

	double Bcrit = 4.4e13;
	double b = 4.0 * gamma * magf/Bcrit;
	//wmag = wmag * fKN_sy(b);


	double deSyn = (4.0/3.0)*thomson*cLight*wmag*	P2(electronMass/p.mass) * gamma*gamma; 
													//P3(electronMass/mass)*(1/(electronMass*cLight2))*(P2(E)/(mass*cLight2));
	return deSyn;
}	


/////////////////////////////////////////////////////////////////
double fSynLosSec(double x, double magneticField, double E)         //funcion a integrar   x=Ega
{
	double mass = electronMass;
	double cte = pow(3.0, 0.5)*P3(electronCharge)*magneticField / (planck*mass*cLight2);

	double Echar = 3 * electronCharge*planck*magneticField*P2(E) / (4 * pi*P3(mass)*cLight*P2(cLight2));

	double tau = x / E;
	double aux = x / (Echar*(1 - tau));  //aca el aux es el x real

	double result = cte*1.85*(1 - tau)*pow(aux, (1.0 / 3.0))*exp(-aux);

	return tau < 1 ? result : 0.0;
} 
///////////////////synchr losses for secondary pairs

double lossesSynSec(double E, double magneticField, Particle& particle)
{
	double EmaxL = 1.6e-12*pow(10.0,17.0);  
	double EminL = 1.6e-12*pow(10.0,-2.0); 

	double result = RungeKuttaSimple(EminL, EmaxL, [&E, &magneticField](double x){
		return fSynLosSec(x, magneticField, E);
	});  //integra entre Emin y Emax
	
	return result;
}
