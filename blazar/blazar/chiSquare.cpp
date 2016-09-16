#include "chiSquare.h"

#include "modelParameters.h"

#include <boost/property_tree/ptree.hpp>

//#include <fparameters\Dimension.h>
//#include <fparameters\SpaceIterator.h>
#include <fparameters/parameters.h>

#include <fmath\physics.h>


double chiSquareFit(ParamSpaceValues psv, Matrix data, int dof)
{
	static const int N = GlobalConfig.get<int>("height");
	//int N = 20;// data.size.size(); //ver cual de las dos es la dimension del vector
	double chi = 0.0;

	for (size_t i = 0; i < N; ++i)
	{
		double aux = 0.0;
		double x_i = data[i][0];
		double observed = data[i][1];
		double err = data[i][2];
		double modelled = psv.interpolate({ { DIM_E, x_i } });
		
		if (err > 0.0) {
			chi = chi + P2((observed - modelled) / err);
		}
	}

	return chi / dof;
	
}

//chequear que las energias del archivo esten dentro de los valores de mi psv, sino se va de rango al interpolar
