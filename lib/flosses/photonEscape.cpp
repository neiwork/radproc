#include "photonEscape.h"

#include <flosses\crossSectionInel.h>
//#include <fparameters\parameters.h>
//#include <fparameters\ParamSpaceValues.h>
#include <fparameters\SpaceIterator.h>
#include <fmath\physics.h>

double photonEscape_f(double x)
{
	if (0.1 < x && x < 1.0){ return (1.0 - x) / 0.9; }
	else if (x <= 0.1){ return 1.0; }
	else { return 0.0; }
}

double escapePhoton(double Eph, Particle& electron, double radius)   //en [s]
{

	double sum = 0.0;
	double dE;
//esto es un promedio de sigma*Ne
	electron.ps.iterate([&dE,&sum,&electron,Eph](const SpaceIterator& i){
		const DimensionIterator& eit = i.its[0]; // energy iterator;
		if (eit.canPeek(1)) {
			double E = i.val(0);//
			dE = eit.peek(0) - i.its[0].peek(1);
			sum = sum + crossSectionKN(Eph, E)*electron.distribution.get(i) * dE;  //E=i.par.E
		}
	});

	//for (size_t i=0; i < electron.energyPoints.size()-1; ++i)
	//{   
	//	dE = electron.energyPoints[i]-electron.energyPoints[i+1];
	//	sum = sum + crossSectionKN(Eph,electron.energyPoints[i])*electron.distribution.values[i]*dE;
	//}

	double tau = sum*radius;

	double Tesc = (radius/cLight)*(1+sum*photonEscape_f(Eph/(electron.mass*cLight2)));

	double original = 2.0*radius/cLight/3.0;

	return Tesc; //original;

	//return 1.0;

}
