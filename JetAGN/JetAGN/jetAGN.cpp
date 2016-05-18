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
#include <fparameters\Dimension.h>
#include <fmath\physics.h>



int jetAGN()
{
	try
	{
		//testBinarySearch();

		setParameters();

		State model;

	//	model.photon.injection.ps.iterate([&model](const SpaceIterator& i){
		//	double nph = blackBody(i.val(DIM_E), i.val(DIM_R));
			//model.nph.set(i, nph);
		//});


	//	radiativeLosses(model);

		injection(model.electron, model);

		checkInyectedPower(model.electron.injection, 0);

	//	writeAllSpaceParam("electronInj.txt", model.electron.injection);
	//	writeEnergyFunction("electronInj_E.txt", model.electron.injection, 1, 1); //escribe Q(E), para r(0) y t(0)
	
		distribution(model.electron, model);
	
		//writeAllSpaceParam("electronDist.txt", model.electron.distribution);
		//writeEandTParamSpace("electronDist_ET.txt", model.electron.distribution, nR/2);
		writeRandTParamSpace("electronDist_RT_2.txt", model.electron.distribution, model.electron.ps[0].size() - 5);
		//writeEnergyFunction("electronDist_E.txt", model.electron.distribution, 1, nR);


		//lo siguiente es una funcion rapida para llenar N(E) asi pruebo las luminosidades
		//model.electron.distribution.fill([&model](const SpaceIterator& i){
		//	return model.electron.injection.get(i);
		//});

		processes(model, "ntLuminosity_LR2.txt");

	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
		throw;
	}

	return 0;
}

