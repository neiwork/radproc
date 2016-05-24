#include "cilindricIntegral.h"

#include "modelParameters.h"
#include <fmath\RungeKutta.h>
#include <fparameters\parameters.h>



double intCilindric(double zMin, double zMax, fun1 fun)
{
	static const double openingAngle = GCFG.get<double>("openingAngle", 0.1);
	const double theta = openingAngle;

	double integral =
		RungeKuttaSimple(zMin, zMax, [theta,fun](double z){
		return fun(z)*(pi*P2(jetRadius(z, theta)));
	});

	return integral;
}




//		RungeKuttaSimple(zMin, zMax,
	///	[theta](double z){
//		return fun(z)*(pi*P2(jetRadius(z, theta)));
//	});
