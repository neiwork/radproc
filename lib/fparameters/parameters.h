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

	int nPhotonEnergies{ 0 };
	double photonLogEmin{ 0 };
	double photonLogEmax{ 0 };	

	// Data of electrons and protons

	double primaryIndex{ 0 };

	double electronLogEmin{ 0 };
	double electronLogEmax{ 0 };

	double protonLogEmin{ 0 };
	double protonLogEmax{ 0 };

	double pionLogEmin{ 0 };
	double pionLogEmax{ 0 };

	double muonLogEmin{ 0 };  //en realidad es log(2*mc^2) = 8.32
	double muonLogEmax{ 0 };

	double neutrinoLogEmin{ 0 };
	double neutrinoLogEmax{ 0 };

	double neutronLogEmin{ 0 };
	double neutronLogEmax{ 0 };

	//Energy and time data

	double E{ 0 };
	double R{ 0 };
	double T{ 0 };
	int nEnergies{ 0 };        //massive particles

	int nTimes{ 0 };
	double timeMin{ 0 };
	double timeMax{ 0 };   //segundos en dos horas

	double rmin{ 0 };
	double rmax{ 0 };
	int nR{ 0 };

};

extern Parameters parameters;