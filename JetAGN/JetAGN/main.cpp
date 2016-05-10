#include <stdio.h>

#include "state.h"
#include "modelParameters.h"
#include "radiativeLosses.h"
#include "injection.h"
////#include "photonInjection.h"
#include "distribution.h"
#include "processes.h"
////#include "photonDistribution.h"

#include "checkPower.h"

#include "write.h"
#include "checks.h"


#include <fparticle\Particle.h>
#include <fparameters\parameters.h>
#include <fmath\physics.h>



int main ()
{
	//testBinarySearch();

	setParameters();

	State model;

//	model.photon.injection.ps.iterate([&model](const SpaceIterator& i){
	//	double nph = blackBody(i.par.E, i.par.R);
		//model.nph.set(i, nph);
	//});


//	radiativeLosses(model);

	injection(model.electron, model);

	checkInyectedPower(model.electron.injection, 0);

//	writeAllSpaceParam("electronInj.txt", model.electron.injection);
//	writeEnergyFunction("electronInj_E.txt", model.electron.injection, 1, 1); //escribe Q(E), para r(0) y t(0)
	
	distribution(model.electron, model);

	
	//writeAllSpaceParam("electronDist.txt", model.electron.distribution);
	writeEandTParamSpace("electronDist_ET.txt", model.electron.distribution, nR/2);
	writeRandTParamSpace("electronDist_RT.txt", model.electron.distribution, nEnergies/2);
	//writeEnergyFunction("electronDist_E.txt", model.electron.distribution, 1, nR);


	//lo siguiente es una funcion rapida para llenar N(E) asi pruebo las luminosidades
	//model.electron.distribution.fill([&model](const SpaceIterator& i){
	//	return model.electron.injection.get(i);
	//});

	//processes(model, "ntLuminosity_LR2.txt");

	return 0;
}

