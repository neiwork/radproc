#include "cilindricIntegral.h"

#include "nonThermalLuminosity.h"

#include "write.h"

#include "modelParameters.h"
#include <fmath\RungeKutta.h>
#include <fparameters\parameters.h>

#include <boost/property_tree/ptree.hpp>


double intCilindric(double zMin, double zMax, fun1 fun)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	const double theta = openingAngle;

	int n = 500;

	double z_int = pow((zMax / zMin), (1.0 / n));

	double z = zMin;
	
	//std::ofstream file;
	//file.open("Mdot_pseda.txt", std::ios::out);
	

	double L1 = 0.0;

	for (int i = 0; i < n; ++i)
	{
		double dz = z*(z_int - 1.0);

		L1 = L1 + fun(z)*(pi*P2(jetRadius(z, theta)))*dz; 

		//file << z / pc << '\t' << L1 << std::endl;
		//file << z / pc << '\t' << L1*yr/solarMass << std::endl;

		z = z*z_int;

	}

	//file.close();

	return L1;

}


//	double integral =
//		RungeKuttaSimple(zMin, zMax, [&](double z){
//		return fun(z)*(pi*P2(jetRadius(z, theta)));
//	});
//
//	return integral;
//}




//		RungeKuttaSimple(zMin, zMax,
	///	[theta](double z){
//		return fun(z)*(pi*P2(jetRadius(z, theta)));
//	});
