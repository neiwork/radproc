#pragma once

#include <fmath\mathematics.h>

class Parameters {
public:

	// <derivados>

		double magneticField{ 0 };
		
		double radius{ 0 };  // cambie Rc por radius, asi que no van a compilar las cosas viejas
		                // [av] usado solo por opticalDepthSSA
	// </derivados>

	double accEfficiency{ 0 };
	double starT{ 0 };

	//jet data
	double openingAngle{ 0 };
	double Gamma{ 0 };
	double Dlorentz{ 0 };
	double Lj{ 0 };

	// power law config
	double primaryIndex{ 0 };

};

extern Parameters parameters;