#include "annihilationLine.h"


#include <iostream>
#include <fstream>
#include <string>

//#include "photonInjection.h"
//#include "photonEscape.h"

#include <finjection\pairAnnihilation.h>
#include <finjection\pairInjection.h>
#include <fparameters\parameters.h>
#include <fmath\physics.h>
#include <iostream>
using namespace std;

void annihilationLine(Particle& photon, Particle& electron, Particle& secondaryElectron, Particle& positron) 
{

	std::ofstream arch;  arch.open("annihilationLine.txt");                                                       

	//for (size_t i=0; i < photon.energyPoints.size(); ++i)
	//{   
	//	double factor = P2(photon.energyPoints[i])*volume;   //este factor pasa de inyeccion a luminosidad a 

	//	double line = pairAnnihilation(photon.energyPoints[i],electron,secondaryElectron,positron);

	//	arch << log10(photon.energyPoints[i]/1.6e-12);

	//	arch << "\t" <<log10(line*factor) << std::endl; 
	//}

	// nueva version (a medio terminar)
	//arch.open("annihilationLine2.txt");

	photon.ps.iterate([&photon,&electron,&secondaryElectron,&positron,&arch](const SpaceIterator& i){

		double factor = P2(i.val(DIM_E))*volume;   //este factor pasa de inyeccion a luminosidad a 

		double line = pairAnnihilation(i.val(DIM_E), electron, secondaryElectron, positron);

		arch << log10(i.val(DIM_E) / 1.6e-12);

		arch << "\t" << log10(line*factor) << std::endl;

	});

}