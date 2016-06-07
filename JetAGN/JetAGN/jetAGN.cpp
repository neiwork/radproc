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
#include "ioutil.h"


#include <fparticle\Particle.h>
#include <fparameters\parameters.h>
#include <fparameters\Dimension.h>
#include <fmath\physics.h>

int jetAGN()
{
	std::string folder{ prepareOutputfolder() };

	try {
		boost::property_tree::ptree cfg{ readConfig() };

		setParameters(cfg);
		State model(cfg);

	//	model.photon.injection.ps.iterate([&model](const SpaceIterator& i){
		//	double nph = blackBody(i.val(DIM_E), i.val(DIM_R));
			//model.nph.set(i, nph);
		//});


	//	radiativeLosses(model);

		injection(model.electron, model);

		std::cout << "checking injected power" << '\t' << computeInjectedPower(model.electron.injection, 0) << std::endl;		

	//	writeAllSpaceParam("electronInj.txt", model.electron.injection);
	//	writeEnergyFunction("electronInj_E.txt", model.electron.injection, 1, 1); //escribe Q(E), para r(0) y t(0)
	
		distribution(model.electron, model);
	
		writeAllSpaceParam(getFileName(cfg, folder, "electronDist"), model.electron.distribution);
		writeEandTParamSpace(getFileName(cfg, folder, "electronDist_ET"), model.electron.distribution, model.electron.ps[1].size() / 2);
		writeRandTParamSpace(getFileName(cfg,folder,"electronDist_RT"), model.electron.distribution, model.electron.ps[0].size()/2);
		writeEnt(getFileName(cfg, folder, "E_NT_r"), model.electron.distribution);


		//lo siguiente es una funcion rapida para llenar N(E) asi pruebo las luminosidades
		//model.electron.distribution.fill([&model](const SpaceIterator& i){
		//	return model.electron.injection.get(i);
		//});

		processes(model, getFileName(cfg,folder,"luminosity"));

	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
		throw;
	}

	return 0;
}

