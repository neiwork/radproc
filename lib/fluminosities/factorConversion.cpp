#include "factorConversion.h"

#include <fparameters\parameters.h>
#include <fmath\physics.h>

double factorLumToNph(double E)
{
	return radius/(volume*cLight*P2(E));  //la diferencia con el otro esta juto en el crossing time
}


double factorLumToQph(double E)
{
	return 1.0/(volume*P2(E));
}