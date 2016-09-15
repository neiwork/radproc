#include "chiSquare.h"

#include "modelParameters.h"

#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters/parameters.h>

#include <fmath\physics.h>

double chiSquare(const ParamSpaceValues psv, Matrix data, int dof)
{
	int N = data.size.size(); //ver cual de las dos es la dimension del vector
	double chi = 0.0;

	for (size_t i = 0; i < N; ++i)
	{
		double x_i = data[i][0];
		double observed = data[i][1];
		double err = data[i][2];
		double modelled = psv.interpolate({ { DIM_E, x_i } });

		chi = chi + P2((observed - modelled) / err);
	}

	return chi / dof;






}
