#include "luminosityHadronic.h"


#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\interpolation.h>
#include <flosses\crossSectionInel.h>
#include <fmath\physics.h>
#include <algorithm>

double fHadron(double x, const Particle& creator,
	const ParamSpaceValues& denf, const SpaceCoord& psc) //funcion a integrar   x=Ecreator; L=L(Ega)
{	

	double Kpi = 0.17;

	double eval = creator.mass*cLight2+x/Kpi;

	double Ekin = x/Kpi;

	double distCreator = creator.distribution.interpolate({ { 0, x } }, &psc);


	double thr = 0.0016; //1GeV

	double sigma = 30e-27*(0.95+0.06*log(Ekin/thr));

	const double density = denf.get(psc);

	double pionEmiss = cLight*density*sigma*distCreator/Kpi;  //sigma = crossSectionHadronicDelta(Ekin)
															  //lo saco asi pongo la condicion Ekin > Ethr en el limite de la int

	double result = pionEmiss/(pow(P2(x)-P2(chargedPionMass*cLight2),0.5));

	return result;
}


double luminosityHadronic(double E, const Particle& creator,
	const ParamSpaceValues& denf, const SpaceCoord& psc)
{
	double Kpi = 0.17;
	double thr = 0.0016; //1GeV

	double Max  = 1.6e-12*pow(10.0,17.0);   //esto es un infinito 
	double Min  = std::max(E+P2(chargedPionMass*cLight2)/(4*E),thr*Kpi); //== Ekin > Ethr

	double integral = RungeKuttaSimple(Min, Max, 
		[&](double x) {return fHadron(x, creator, denf, psc); }
	);    //integra entre Emin y Emax

	double luminosity = integral*P2(E); // [erg s^-1 cm^-3 ]

	return luminosity; //P2(Dlorentz);  

	//Dlorentz es el factor que transforma las distribuciones en el caso de jets
	// con /P2(Dlorentz) paso del sist de lab al sist comoving con el jet ;   
}

