#include "ppPionInj.h"

#include "dataInjection.h"
#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <flosses\crossSectionInel.h>
#include <fmath\interpolation.h>
	

double fPP(double x, double E, Particle& creator )         //funcion a integrar   x=Eproton; E=Epion
{
	//DataInjection* data = (DataInjection*)voiddata;
	//double E = data->E;    //esta E corresponde a la energia del foton emitido; E=Ega
	//double mass   = data->mass;
	//Vector& Ncreator = data->Ncreator;
	//Vector& Ecreator = data->Ecreator;

	//const double mass   = particle.mass;
	//const Vector& Ncreator = creator.distribution.values;
	//const Vector& Ecreator = creator.eDim()->values;


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

	//double distProton = interpol(x,Ecreator,Ncreator,Ncreator.size()-1);
	double distProton = creator.dist(x);

	return f*distProton*se/x;		//Q   = f(Ep,x,ap)*se*Np*dEp/Ep
}

double ppPionInj(double E, Particle& proton)
{                                                                    //el proton es el creator del pion
	double Emax = 1.6e-12*pow(10.0, proton.logEmax);

	//DataInjection data;

	//data.E = E;
	//data.mass = particle.mass;
	//data.Ncreator = proton.distribution.values;
	//data.Ecreator = proton.energyPoints;

	double integralP = RungeKuttaSimple(E, Emax, [&E, &proton](double x){
		return fPP(x, E, proton);
	});

	double injection = integralP*cLight*density;

	return injection;   ///P2(Dlorentz);  
	//Dlorentz es el factor que transforma las distribuciones en el caso de jets

	// con /P2(Dlorentz) paso del sist de lab al sist comoving con el jet
}   //lo comento porque ya lo estoy calculando en el comovil

