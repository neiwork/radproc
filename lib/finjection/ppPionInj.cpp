#include "ppPionInj.h"

#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <flosses/crossSectionInel.h>
//#include <fmath\interpolation.h>
	

double fPP(double Ep, double E, const Particle& creator, const SpaceCoord& psc)         //funcion a integrar   x=Eproton; E=Epion
{ 
	double L       = log(Ep/1.602); //el 1.6 son TeV en erg
	double ap      = 3.67+L*(0.83+L*0.075);
	double sigmapp = crossSectionHadronic(Ep);
	double Bpi     = ap+0.25;
	double r       = 2.6/sqrt(ap);
	double alpha   = 0.98/sqrt(ap);
	double x  	   = E/Ep;
	double xalpha = pow(x,alpha);
	double factor  = 1.0-xalpha;

	double f = 0.0;
	if (factor > 0){
		f = 4.0*alpha*Bpi * xalpha/x * pow(factor/(1.0+r*xalpha*factor),4) * 
				(1.0/factor + r*(1.0-2.0*xalpha)/(1.0+r*xalpha*factor))*sqrt(1.0-chargedPionMass*cLight2/(x*Ep));
		}
	//else
	//	f = 0.0;
	
	double distProton = creator.distribution.interpolate({{0,Ep}},&psc);
	return f*distProton*sigmapp/Ep;
}

//#include <iostream>
//double ppPionInj(double E, Particle& proton
double ppPionInj(double E, const Particle& creator, const double density, const SpaceCoord& psc)
{
	double Emax = creator.emax();
	double integralP = integSimpson(log(E),log(Emax),[E,&creator,&psc](double logEp)
						{
							double Ep = exp(logEp);
							return Ep*fPP(Ep,E,creator,psc);
						},50);

	//double integralP = RungeKuttaSimple(E, Emax, [&E, &creator, &psc](double x){
	//	return fPP(x, E, creator, psc);
	//});

	double injection = integralP*cLight*density;

	return injection;  
}


