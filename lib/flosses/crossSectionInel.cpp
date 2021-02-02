#include "crossSectionInel.h"


#include <fparticle/Particle.h>
#include <fparameters/parameters.h>
#include <fmath/physics.h>

double crossSectionPHPion(double t)
{
	double aux1 = 500.0e6*1.602e-12;  //500MeV

	if (pionThresholdPH <= t && t <= aux1)	{
		return 340.0e-30;  //microb
	}
	else if ( t > aux1)	{
		return 120.0e-30;  //microb
	}
	else	{
		return 0;
	}
}

double inelasticityPHPion(double t)
{	
	double aux1 = 500.0e6*1.602e-12;  //500MeV

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
	double aux2	= 4.0*electronRestEnergy;
	double x = t / electronRestEnergy;

	if (pairThresholdPH <= t && t <= aux2) {

		//return 1.2135e-27*pow((x-2)/2,3);  //cm2
		double eta = (x-2.0)/(x+2.0);
		return 0.25*thomson*fineStructConst*pow((x-2)/x,3)*(1.0+0.5*eta+23.0/40.0*eta*eta+
					37.0/120.0*eta*eta*eta+61.0/192.0*eta*eta*eta*eta);  //cm2
	} else if ( t > aux2) {
		return fineStructConst*electronRadius*electronRadius*(3.1111*log(2*x)-8.0741+
			   P2(2/x)*(6.0*log(2*x)-3.5+2.0/3.0*P3(log(2*x))-P2(log(2*x))-1.0/3.0*pi*pi*log(2*x)
						+pi*pi/6.0));

		//return 5.7938e-28*(3.1111*log(2*x)-8.0741+
		//	   P2(2/x)*(2.7101*log(2*x)-P2(log(2*x))+0.6667*pow(log(2*x),3)+0.5490));  //cm2
	}
	else	{
		return 0.0;
	}
}

double inelasticityBetheHeitler(double t)
{	
	double aux3	=	1000.0;
	double x = t/electronRestEnergy;

	double factor = 4.0*(electronMass/protonMass)/x;

	if (pairThresholdPH <= t && x <= aux3)	{
		return factor*(1+0.3957*log(x-1)+0.1*P2(log(x-1))+0.0078*P3(log(x-1)));  
	}
	else if ( x > aux3)	{
		return factor*((-8.78+5.513*log(x)-1.612*P2(log(x))+0.668*P3(log(x)))/
			           (3.1111*log(2*x)-8.0741));
	}
	else	{
		return 0.0;
	}
}

double crossSectionHadronic(double E)
{
	//double L = log(E/1.6);  //el 1.6 son TeV en erg
	double Tp = E-protonMass*cLight2;
	double Tp_th = pionThresholdH - protonMass*cLight2;
	//double sigmaHadronic = 1.0e-27*(34.3+1.88*L+0.25*P2(L))*P2(1-pow((pionThresholdH/E),4));  //mb	
	double sigmaHadronic = 1.0e-27 * (30.7-0.96*log(Tp/Tp_th)+0.18*P2(log(Tp/Tp_th))) * 
							P3(1.0-pow(Tp_th/Tp,1.9));
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


double angleAveragedKN(double eps)
{
	double l = log(1.0+2.0*eps);
	double sigmaKN = (3.0*thomson/4.0) *
					( (1.0+eps)* (2.0*eps*(1.0+eps)/(1.0+2.0*eps) - l) / P3(eps) 
	               + l/(2.0*eps) - (1.0+3.0*eps)/P2(1.0+2.0*eps) );
		
	double result = 0.0;
	
	if(eps < 1)
	{
		result = thomson*(1.0-2.0*eps+26.0*eps*eps/5.0);
	}
	else
	{
		result = 3.0*thomson*(log(2.0*eps)+0.5)/(8.0*eps);
	}
		
	return result;
}

double crossSectionGammaGamma(double beta)
{
	return  3.0*thomson*(1.0-P2(beta))*
		    (2.0*beta*(P2(beta)-2.0)+
			(3.0-pow(beta,4))*log((1.0+beta)/(1.0-beta)));
}