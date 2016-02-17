#include "crossSectionInel.h"

#include <fmath\RungeKutta.h>
#include <fparticle\particle.h>
#include <fparameters\parameters.h>
#include <fmath\physics.h>

double crossSectionPHPion(double t)
{
	double aux1 = 500.0e6*1.6e-12;  //500MeV

	if (pionThresholdPH <= t && t <= aux1)	{
		return 340.0e-30;  //mb
	}
	else if ( t > aux1)	{
		return 120.0e-30;  //mb
	}
	else	{
		return 0;
	}
}

double inelasticityPHPion(double t)
{	
	double aux1 = 500.0e6*1.6e-12;  //500MeV

	if (pionThresholdPH <= t && t <= aux1)	{
		return 0.2;  
	}
	else if ( t > aux1)	{
		return 0.6;  
	}
	else	{
		return 0.0;
	}
}

double crossSectionBetheHeitler(double t)
{
	double aux2	= 4.0*electronMass*cLight2;

	double x	= t/(electronMass*cLight2);

	if (pairThresholdPH <= t && t <= aux2)	{
		return 1.2135e-27*pow((x-2)/2,3);  //cm2
	}
	else if ( t > aux2)	{
		return 5.7938e-28*(3.1111*log(2*x)-8.0741+
			   P2(2/x)*(2.7101*log(2*x)-P2(log(2*x))+0.6667*pow(log(2*x),3)+0.5490));  //cm2
	}
	else	{
		return 0.0;
	}
}

double inelasticityBetheHeitler(double t)
{	
	double mass = protonMass;

	double aux3	=	1000.0*electronMass*cLight2;

	double x = t/(electronMass*cLight2);

	double factor = 4.0*(electronMass/mass)/x;

	if (pairThresholdPH <= t && t <= aux3)	{
		return factor*(1+0.3957*log(x-1)+0.1*P2(log(x-1))+0.0078*pow(log(x-1),3));  
	}
	else if ( t > aux3)	{
		return factor*((-8.78+5.513*log(x)-1.612*P2(log(x))+0.668*pow(log(x),3))/
			           (3.1111*log(2*x)-8.0741));
	}
	else	{
		return 0.0;
	}
}

double crossSectionHadronic(double E)
{
	double L = log(E/1.6);  //el 1.6 son TeV en erg
	
	double sigmaHadronic = 1.0e-27*(34.3+1.88*L+0.25*P2(L))*P2(1-pow((pionThresholdH/E),4));  //mb	

	return E > pionThresholdH ? sigmaHadronic : 0.0;
}

double crossSectionHadronicDelta(double Ekin)
{
	double thr = 0.0016; //1GeV

	double sigma = 30.0e-27*(0.95+0.06*log(Ekin/thr));

	return Ekin > thr ? sigma : 0.0;
	
}

double crossSectionThomson(double mass)
{
	return 8.0*pi*P2(P2(electronCharge)/(mass*cLight2))/3.0;
}

double crossSectionKN(double Eph, double Ee)
{
	double sigmaKN = 3.0*thomson*P2(electronMass*cLight2)*log(2.0*Ee*Eph/P2(electronMass*cLight2)+0.5)/(8.0*Ee*Eph);
	return sigmaKN;
}

double crossSectionGammaGamma(double beta)
{
	return  3.0*thomson*(1.0-P2(beta))*
		    (2.0*beta*(P2(beta)-2.0)+
			(3.0-pow(beta,4))*log((1.0+beta)/(1.0-beta)));
}