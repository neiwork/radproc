#include "nonThermalLuminosity.h"

#include "dynamics.h"
//#include "modelParameters.h"

#include "targetFields.h"
//#include "cilindricIntegral.h"

//#include "write.h"

#include <flosses\nonThermalLosses.h>
//#include <fparameters\SpaceIterator.h>
//#include <fparameters\Dimension.h>
#include <fparameters\parameters.h>


#include <fmath\physics.h>

#include <iostream>
#include <boost/property_tree/ptree.hpp>



double frad(double E, double z, double Gc) 
{
	
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	
	double vel_lat = cLight*openingAngle;

	double tad = adiabaticLosses(E, z, vel_lat, Gc) / E; //FF
			
	double wph = wphIR(z,"IR")*P2(Gc);

	double tic = 4.0*cLight*thomson*wph*(E / P2(electronMass*cLight2)) / 3.0;

	double wmag = P2(computeMagField(z, Gc)) / (8.0*pi);
	double tsyn = 4.0*cLight*thomson*wmag*(E / P2(electronMass*cLight2)) / 3.0;
	
	double frad = tic / (tad + tic + tsyn);
	//double frad = tsyn / (tad + tic + tsyn);
	return frad;

}


double dLnt(double z, double Gc, double z_int, double Rs)
{
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double theta = GlobalConfig.get<double>("openingAngle");
	
	/*double y = z / z_int;

	double g = Gc / Gj;
	
	double frac = (Rs/jetRadius(z, theta));

	if (frac > 1.0) {
		Rs = jetRadius(z, theta); 
	//	std::cout << "SN radius larger than jet radius" << '\t' << 	frac 
	//		<< '\t' << z/pc << std::endl;
	}*/

	//double Lsc = accEfficiency*(1.0/P2(g)-P2(g))*cte/(4.0*P2(Gj));  //mhd
	//double Lsc = accEfficiency*P2(1.0 / g - g)*cte / (4.0*P2(Gj));  //hydro
	
	double beta_j = beta(Gj);
	double beta_c = beta(Gc);
	double cte2 = accEfficiency*jetRamPress(z)*pi*cLight*P2(Rs); // =eta*Lj*So/Sj
	double Lsc = cte2*beta_c*Gc*P2(1.0 - beta_c / beta_j);  //hydro

	return Lsc;
}


