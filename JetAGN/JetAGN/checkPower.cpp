#include "checkPower.h"

#include <fparameters/ParamSpaceValues.h>
#include <fparameters/ParamSpace.h>
#include <fparameters/Dimension.h>

#include "write.h"

#include <fmath\physics.h>


double computeInjectedPower(const ParamSpaceValues& dist, int t_ix)
{
	double sum = 0.0;

	double T = dist.ps[2][t_ix];

	Vector& z = dist.ps.dimensions[1]->values;
	Vector& E = dist.ps.dimensions[0]->values;

	for (size_t i = 0; i < z.size() - 1; ++i) { //no llego al ultimo
		
		double dz = z[i + 1] - z[i];

		double jetR = jetRadius(z[i], parameters.openingAngle);

		for (size_t j = 0; j < E.size() - 1; ++j) {

			double dE = E[j + 1] - E[j];

			double emissivity = dist.interpolate({ { DIM_E, E[j] }, { DIM_R, z[i] }, { DIM_T, T } });

			double L1 = 2.0*pi*P2(jetR)*emissivity*E[j] * dE*dz;

			sum = sum + L1;
		}

	}

	return sum;
}