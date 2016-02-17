#include "luminosityHadronic.h"

#include "packageData.h"
#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\interpolation.h>
#include <flosses\crossSectionInel.h>
#include <fmath\physics.h>
#include <algorithm>

double fHadron(double x, Particle& creator)         //funcion a integrar   x=Ecreator; L=L(Ega)
{
	//Data* data = (Data*)voiddata;
	//double E = data->E;                       //esta E corresponde a la energia del foton emitido; E=Ega
	//double mass = data->mass;
	//Vector& Ncreator = data->Ncreator;
	//Vector& Ecreator = data->Ecreator;	

	double Kpi = 0.17;

	double eval = creator.mass*cLight2+x/Kpi;

	double Ekin = x/Kpi;

	double distCreator = creator.dist(eval);

	double thr = 0.0016; //1GeV

	double sigma = 30e-27*(0.95+0.06*log(Ekin/thr));

	double pionEmiss = cLight*density*sigma*distCreator/Kpi;  //sigma = crossSectionHadronicDelta(Ekin)
															  //lo saco asi pongo la condicion Ekin > Ethr en el limite de la int

	double result = pionEmiss/(pow(P2(x)-P2(chargedPionMass*cLight2),0.5));

	return result;
}


double luminosityHadronic(double E, Particle& creator)
{
	double Kpi = 0.17;
	double thr = 0.0016; //1GeV

	double Max  = 1.6e-12*pow(10.0,17.0);   //esto es un infinito 
	double Min  = std::max(E+P2(chargedPionMass*cLight2)/(4*E),thr*Kpi); //aca pongo la condicion Ekin > Ethr

	//Data data;

	//data.E = E;
	//data.mass = creator.mass;
	//data.Ncreator = Ncreator;
	//data.Ecreator = creator.energyPoints;
 //
	double integral = RungeKuttaSimple(Min, Max, [&creator](double x){
		return fHadron(x,creator);
	});    //integra entre Emin y Emax

	double luminosity = integral*P2(E)*volume; //multiplico por E asi obtengo luminosidad

	return luminosity; //P2(Dlorentz);  

	//Dlorentz es el factor que transforma las distribuciones en el caso de jets
	// con /P2(Dlorentz) paso del sist de lab al sist comoving con el jet ;   
}



/* 	if (x > 0.1e12*1.6e-12){

		double distCreator = interpolMod(x,Ecreator,Ncreator,Ncreator.size()-1);

		double L = log(E/1.6);  //el 1.6 son TeV en erg

		double equis = E/x;

		double Bga  = 1.30+0.14*L+0.011*P2(L);

		double beta = pow(1.79+0.11*L+0.008*P2(L),-1);

		double factor = pow(equis,beta);

		double kga  = pow(0.801+0.049*L+0.014*P2(L),-1); 

		double f = Bga*(log(equis)/equis)*pow((1-factor)/(1+kga*factor*(1-factor)),4)*
				   (1/log(x)-4*beta*factor/(1-factor)-4*kga*beta*factor*(1-2*factor)/(1+kga*factor*(1-factor)));
 
		result = cLight*density*crossSectionHadronic(x)*distCreator*f/x;
	}
	else {*/
