#include "lossesIC.h"

#include "crossSectionInel.h"
//#include "dataLosses.h"
#include <fmath\RungeKutta.h>
#include <fparameters\parameters.h>
#include <fmath\physics.h>

double cIC(double u)   //limite inferior
{
	return u;
}

double dIC(double u, double E, double mass)   //limite superior     
{

	double s = 4*u*E/P2(mass*cLight2);
	return s*E/(1+s);
}

double fIC(double u,double t, double E, double mass, fun1 tpf)   //funcion a integrar
{

	double Erep = mass*cLight2;
	double r    = t*P2(Erep)/(4*u*E*(E-t));
	double Nterm = tpf(u);

	double result = (Nterm/u)*(t-u)*(2*r*log(r)+
       	(1+2*r)*(1-r)+(P2((t/(E-t)))*(1-r))/(2*(1+(t/(E-t)))));

	return  result;
}

double lossesIC(double E, Particle& particle, fun1 tpf, double phEmin, double phEmax)
{
	double mass = particle.mass;

	double constant  = 3*crossSectionThomson(mass)*P2(mass)*pow(cLight,5)/4;

	double a  = phEmin;     //energia minima de los fotones en erg
	double b  = phEmax;     //energia maxima de los fotones en erg

	double integral = RungeKutta(a,b,&cIC,
		[mass, E](double u){return dIC(u, E, mass);},
		[E, mass, tpf](double u, double t){ return fIC(u, t, E, mass, tpf); });    //le asigno a la variable integral el resultado de la integracion   

	double gamma = E/(mass*cLight2);

	double de = constant*integral/P2(E);
			
	return de;

	}
