#include "lossesBrem.h"

#include <fparameters\parameters.h>
#include <fmath\physics.h>

double lossesBremss(double E, double density, Particle& particle)
{
	double mass = particle.mass;
	

	double factorD = 1;
	double Elab = E/factorD;      //E es la energía en el sist comovil=jet

	double deBrem = 4*density*P2(electronRadius)*fineStructConst*cLight*
		            (log(2*Elab/(mass*cLight2))-1.0/3.0)*E/factorD;
	                                                   //divido por factor porque transformo t^-1
	                                                   //multiplico por E y no Elab porque ya transforme t y
                                                       //E es la del comovil 

	return deBrem;
}