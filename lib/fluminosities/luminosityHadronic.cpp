#include "luminosityHadronic.h"


#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\interpolation.h>
#include <flosses\crossSectionInel.h>
#include <fmath\physics.h>
#include <algorithm>

double fHadron(double x, const Particle& p, const double density, const SpaceCoord& psc) //funcion a integrar   x=Ecreator; L=L(Ega)
{	
	double Kpi = 0.17;
	double eval = p.mass*cLight2+x/Kpi;
	
	//double Ekin = Ep/Kpi;
	
	double distCreator=0.0;
	if (eval < p.emax() && eval >= p.emin()) {
		distCreator = p.distribution.interpolate({ { 0, eval } }, &psc);
	}
	
	//double thr = 0.0016; //1GeV
	//double sigma = 30e-27*(0.95+0.06*log(Ekin/thr));
	
	double l = log10((protonMass*cLight2+x/Kpi)/1.6); //evaluada en eval
	double sigma = 1.e-27 * (34.3+1.88*l+0.25*l*l);
	double pionEmiss = cLight*density*sigma*distCreator/Kpi;  //sigma = crossSectionHadronicDelta(Ekin)
															  //lo saco asi pongo la condicion Ekin > Ethr en el limite de la int
	
	double result = 2.0*pionEmiss/(pow(P2(x)-P2(chargedPionMass*cLight2),0.5));
	return result;
}


double luminosityHadronic(double E, const Particle& creator,
	const double density, const SpaceCoord& psc)
{
	double Kpi = 0.17;
	double thr = 0.0016; //1GeV

	double Max  = 1.6e-12*pow(10.0,17.0);   //esto es un infinito 
	double Min  = std::max(E+P2(chargedPionMass*cLight2)/(4*E),thr*Kpi); //== Ekin > Ethr

	double integral = RungeKuttaSimple(Min, Max, 
		[&](double x) {return fHadron(x, creator, density, psc); }
	);    //integra entre Emin y Emax

	double luminosity = integral*E*E; // [erg s^-1 cm^-3 ]

	return luminosity; 
}

