#include "ppPionInj.h"

#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <flosses\crossSectionInel.h>
#include <fmath\interpolation.h>
	

double fPP(double x, double E, const Particle& creator, const SpaceCoord& psc)         //funcion a integrar   x=Eproton; E=Epion
{ 
	double L      = log(x/1.6); //el 1.6 son TeV en erg
	double ap     = 3.67+0.83*L+0.075*P2(L);
	double se     = crossSectionHadronic(x);
	double Bp     = ap+0.25;
	double r      = 2.6*pow(ap,-0.5);
	double alpha  = 0.98*pow(ap,-0.5);
	double equis  = E/x;
	double factor = 1-pow(equis,alpha);

	double f;
	if (factor =! 0)	{
		f      = 4*alpha*Bp*pow(equis,(alpha-1))*pow((factor/(1+r*pow(equis,alpha)*factor)),4)
      	         *(1/factor+r*(1-2*pow(equis,alpha))/(1+r*pow(equis,alpha)*factor))*
     	         pow((1-chargedPionMass*cLight2/(equis*x)),0.5);
	}
	else	{
		f = 0.0;
	}

	double distProton = creator.distribution.interpolate({ { 0, x } }, &psc); // creator.dist(x);

	return f*distProton*se/x;		//Q   = f(Ep,x,ap)*se*Np*dEp/Ep
}

//double ppPionInj(double E, Particle& proton
double ppPionInj(double E, const Particle& creator,
		const double density, const SpaceCoord& psc)
{                                                                    //el proton es el creator del pion
	double Emax = 1.6e-12*pow(10.0, creator.logEmax);


	double integralP = RungeKuttaSimple(E, Emax, [&E, &creator, &psc](double x){
		return fPP(x, E, creator, psc);
	});

	double injection = integralP*cLight*density;

	return injection;   
}

