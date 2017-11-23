#include "targetFields.h"



#include "modelParameters.h"
#include "functions.h"
#include <fparameters/SpaceIterator.h>

//#include <fmath/physics.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>


void magFieldFill(State& st)
{
	static const double beta = GlobalConfig.get<double>("beta");

	st.magf.fill([&](const SpaceIterator& i) {
		
		double r = i.val(DIM_R);
	    double theta = i.val(DIM_THETA);
		return sqrt( beta * 24.0 * pi * pressureTot(r, theta) );
	});
}

	  

