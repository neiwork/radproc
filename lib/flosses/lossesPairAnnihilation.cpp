#include "lossesPairAnnihilation.h"

#include "dataLosses.h"
#include <fparameters\parameters.h>
#include <flosses\crossSectionInel.h>
#include <fmath\interpolation.h>
#include <fmath\RungeKutta.h>
#include <fmath\physics.h>

double fAnn(double x, double E, Particle& secondaryParticle)         //funcion a integrar   
{
	//DataLosses* data = (DataLosses*)voiddata;
	//double E = data->E;                    //esta E corresponde a la particula que esta perdiendo energia
	//Vector& Ntarget(data->Ntarget);
	//Vector& Etarget(data->Etarget);

	double Erest = electronMass*cLight2;
	double cte = 3.0*thomson*cLight*P2(Erest) / 8.0;

	double distTarget = secondaryParticle.dist(x);
	//interpol(x, Etarget, Ntarget, Etarget.size() - 1);

	return cte*distTarget*(log(4 * x*E / P2(Erest)) - 2) / (x*E);
}

double lossesPairAnnihilation(double E, Particle& particle, Particle& secondaryParticle)
{
	double Emax = 1.6e-12*pow(10.0, particle.logEmax);
	double Emin = 1.6e-12*pow(10.0, particle.logEmin);

	double integral = E*RungeKuttaSimple(Emin, Emax, [E,&secondaryParticle](double x){
		return fAnn(x, E, secondaryParticle);
	});    //integra entre Emin y Emax
	return integral;                                                //*E para que sea dE/dt
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/*double fAnnEl(double x, void* voiddata)         //funcion a integrar   
{
	lossesData* data = (lossesData*)voiddata;
	double E = data->E;                    //esta E corresponde a la particula que esta perdiendo energia
	Vector& Np(data->current->Np);
	Vector& Ep(data->current->Ep);  //????

	double Erest = electronMass*cLight2;
	double cte = 3*thomson*cLight*P2(Erest)/8.0;
	
	double distPositron = interpol(x,Ep,Np,Ep.size()-1);

	return  cte*distPositron*(log(4*x*E/P2(Erest))-2)/(x*E);  
}

double LossesAnnEl(double E, const CoupledEqSys* previous, CoupledEqSys* current)
{
	double Emax = 1.6e-12*pow(10.0,protonLogEmax);  
	double Emin = 1.6e-12*pow(10.0,electronLogEmin); 

	lossesData ldata;

	ldata.E = E;
	ldata.current = current;
	ldata.previous = previous;

	double integral = E*RungeKuttaSimple(Emin, Emax, fAnnEl, &ldata);    //integra entre Emin y Emax
	return integral ;   
}*/ 
