#include "checkPower.h"

#include "write.h"

#include <fparameters\SpaceIterator.h>

#include <fparameters/ParamSpaceValues.h>
#include <fparameters/ParamSpace.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>

#include <fmath\physics.h>

#include <boost/property_tree/ptree.hpp>




double computeInjectedPower(const ParamSpaceValues& dist)
{
	
	
	const double EMIN = dist.ps[DIM_E].first();
	const double EMAX = dist.ps[DIM_E].last();
	const int N_E = dist.ps[DIM_E].size() - 1;
	double E_int = pow((EMAX / EMIN), (1.0 / N_E));

	//const double RMIN = dist.ps[DIM_R].first();
	//const double RMAX = dist.ps[DIM_R].last();
	//const int N_R = dist.ps[DIM_R].size() - 1;
	//double sumT = 0.0;
	
	//for (int z_ix = 0; z_ix < N_R; z_ix++) {
		
		double sum = 0.0;

		double E = EMIN;

		//double z = dist.ps[DIM_R][z_ix];

		dist.ps.iterate([&](const SpaceIterator& i) {
			
			double dE = E*(E_int - 1);

			double emissivity = dist.interpolate({ { DIM_E, E }});

			double L1 = emissivity* E * dE;

			sum = sum + L1;

			E = E*E_int;

		});


		std::cout << "checking injected power:" << '\t' << sum << std::endl;

		return sum;

	
}

/*

double computeInjectedEnergy(const ParamSpaceValues& dist)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	
	Vector& E = dist.ps.dimensions[0]->values;

	//for (size_t i = 0; i < z.size() - 1; ++i) { //no llego al ultimo

	//	double dz = z[i + 1] - z[i];

	//	double jetR = jetRadius(z[i], openingAngle);

		for (size_t j = 0; j < E.size() - 1; ++j) {


			double dE = E[j + 1] - E[j];

			double emissivity = dist.interpolate({ { DIM_E, E[j] } });

			//double L1 = 2.0*pi*P2(jetR)*emissivity*E[j] * dE*dz;
			double L1 = emissivity*E[j] * dE;

			//sum = sum + L1;
		}

		
	//}

	return sum;
}*/ 

//Vector& z = dist.ps.dimensions[1]->values;
//Vector& E = dist.ps.dimensions[0]->values;

//	double Emin = p.emin();