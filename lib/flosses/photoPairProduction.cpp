#include "photoPairProduction.h"

#include "crossSectionInel.h"
#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\interpolation.h>
#include <fmath\physics.h>
#include <algorithm>

double fSink(double e2, double E, fun1 tpf)         //funcion a integrar
{

	double Erest = electronMass*cLight2;
	
	double e1   = E/Erest;  //esta es la misma e1 de afuera 
	
	double beta = sqrt(1.0-1.0/(e1*e2));

	double photonDist_y = tpf(e2*Erest);

	double result = photonDist_y*crossSectionGammaGamma(beta);  

	return result;    
}


double photoPairProduction(double E, const Particle& photon, fun1 tpf, double tpEmin, double tpEmax)
{
	//
	double Erest = electronMass*cLight2;
	
	double cte = cLight*Erest;  // *Erest^2: de las dos dist de fotones
	                            // /Erest: de la inyeccion de electrones

	//normalizo todo a me*c^2

	double e1   = E/Erest;

	double uno  = 2.0/e1;     //esta es la condicion sobre las energias de los
	                          //fotones para que puedan crear pares  
	double dos  = tpEmin/Erest;

	double inf = std::max(uno,dos);

	double sup = tpEmax/Erest;  											

	if(inf < sup){
		double integral = RungeKuttaSimple(inf, sup, [E,tpf](double e){
			return fSink(e, E, tpf);
		});

		double photonDist_x = tpf(e1*Erest); 

		double emissivityA = cte*photonDist_x*integral;

		return emissivityA;
	}
	else{ return 0.0;}

}