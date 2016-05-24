#pragma once

#include <fmath\mathematics.h>

class Parameters {
public:

	// <derivados>

		//double magneticField{ 0 };
		
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