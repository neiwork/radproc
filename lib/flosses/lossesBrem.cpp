#include "lossesBrem.h"

#include <fmath/physics.h>

double lossesBremss(double E, double density, Particle& particle)
{
	double mass = particle.mass;
	

	double deBrem = 4*density*P2(electronRadius)*fineStructConst*cLight*
		            (log(2*E/(mass*cLight2))-1.0/3.0)*E;
	                                                  
	return deBrem;
}