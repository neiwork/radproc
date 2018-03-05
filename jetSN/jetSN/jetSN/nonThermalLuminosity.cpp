#include "nonThermalLuminosity.h"

#include "dynamics.h"
//#include "modelParameters.h"

#include "targetFields.h"
//#include "cilindricIntegral.h"

//#include "write.h"

#include "lossesAnisotropicIC.h"
#include <flosses\lossesSyn.h>
#include <flosses\nonThermalLosses.h>
//#include <fparameters\SpaceIterator.h>
//#include <fparameters\Dimension.h>
#include <fparameters\parameters.h>


#include <fmath\physics.h>

#include <iostream>
#include <boost/property_tree/ptree.hpp>




double frad_2(double E, double z, double Gc, Particle& electron)
{

	static const double starTIR = GlobalConfig.get<double>("IRstarT");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	double vel_lat = cLight*openingAngle;

	double tad = adiabaticLosses(E, z, vel_lat, Gc) / E; //FF

	double wph = wphIR(z, "IR")*P2(Gc);

	double tic = 4.0*cLight*thomson*wph*(E / P2(electronMass*cLight2)) / 3.0;

	double B = computeMagField(z, Gc);
	
	double eIC_Aux = lossesAnisotropicIC(E, electron, z, Gc,
		[&](double E, double z, double r) {
		return nph_ICani(E, z, r, Gc, "IR"); },
		starTIR) / E;

		double eSyn = lossesSyn(E, B, electron) / E;

		double frad = eIC_Aux / (tad + eIC_Aux + eSyn);

		//double frad = tic / (tad + tic + tsyn);
		//double frad = tsyn / (tad + tic + tsyn);
		return frad;

}


double frad(double E, double z, double Gc) 
{
	
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	
	double vel_lat = cLight*openingAngle;

	double tad = adiabaticLosses(E, z, vel_lat, Gc) / E; //FF
			
	double wph = wphIR(z,"IR")*P2(Gc);

	double tic = 4.0*cLight*thomson*wph*(E / P2(electronMass*cLight2)) / 3.0;

	double B = computeMagField(z, Gc);
	double wmag = P2(B) / (8.0*pi);
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
	
	double hj = 1.0;
		
	double beta_j = beta(Gj);
	double beta_c = beta(Gc);

	double beta_rel = (beta_j - beta_c) / (1.0 - beta_j*beta_c);
	double G_rel = 1.0 / sqrt(1.0 - P2(beta_rel));


	double cte2 = accEfficiency*jetRamPress(z)*pi*cLight*P2(Rs); // =eta*Lj*So/Sj
	//double Lsc = cte2*beta_c*Gc*P2(1.0 - beta_c / beta_j);  //hydro
	double Lsc = cte2*beta_rel*G_rel*(G_rel*hj - 1.0)/(P2(Gj)*beta_j);  //hydro

	return Lsc;
}


