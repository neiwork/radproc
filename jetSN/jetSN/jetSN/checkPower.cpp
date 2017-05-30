#include "checkPower.h"

#include "write.h"

#include <fparameters\SpaceIterator.h>

#include <fparameters/ParamSpaceValues.h>
#include <fparameters/ParamSpace.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>

#include <fmath\physics.h>

#include <boost/property_tree/ptree.hpp>




double computeInjectedPower(const ParamSpaceValues& dist, int t_ix)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	double sum = 0.0;

	//double T = dist.ps[2][t_ix];

	///Vector& z = dist.ps.dimensions[1]->values;
	Vector& E = dist.ps.dimensions[0]->values;


	//const double RMIN = dist.ps[DIM_R].first();
	//const double RMAX = dist.ps[DIM_R].last();
	//const int N_R = dist.ps[DIM_R].size() - 1;

	//double z_int = pow((RMAX / RMIN), (1.0 / N_R));
	
	//for (size_t i = 0; i < N_R + 1 ; ++i) { 
	
		//double z = z[i]; // i.val(DIM_R);
		//double dz = z[i]*(z_int - 1);
		//volumen de la celda i
		//double vol_i = pi*P2(jetRadius(z[i], openingAngle))*dz;

		//double jetR = jetRadius(z[i], openingAngle);

		for (size_t j = 0; j < E.size() - 1; ++j) {

			double dE = E[j + 1] - E[j];

			double emissivity = dist.interpolate({ { DIM_E, E[j] } });

			//double L1 = 2.0*pi*P2(jetR)*emissivity*E[j] * dE*dz;
			double L1 = emissivity*E[j] * dE;

			sum = sum + L1;
		}

	//}

	return sum;
}


double computeInjectedEnergy(const ParamSpaceValues& dist)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	double sum = 0.0;

	//double T = dist.ps[2].last();

	//Vector& z = dist.ps.dimensions[1]->values;
	Vector& E = dist.ps.dimensions[0]->values;

	//for (size_t i = 0; i < z.size() - 1; ++i) { //no llego al ultimo

	//	double dz = z[i + 1] - z[i];

	//	double jetR = jetRadius(z[i], openingAngle);

		for (size_t j = 0; j < E.size() - 1; ++j) {

			double dE = E[j + 1] - E[j];

			double emissivity = dist.interpolate({ { DIM_E, E[j] } });

			//double L1 = 2.0*pi*P2(jetR)*emissivity*E[j] * dE*dz;
			double L1 = emissivity*E[j] * dE;

			sum = sum + L1;
		}

	//}

	return sum;
}

