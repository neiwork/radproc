#include <stdio.h>


#include "maximumEnergy.h"
#include "checkPower.h"
#include "processes.h"
#include "radiativeLosses.h"
#include "nonThermalLuminosity.h"
#include "distribution.h"
#include "injection.h"
#include "state.h"
#include "modelParameters.h"

#include "write.h"
#include "checks.h"
#include "ioutil.h"


#include <fparticle\Particle.h>
#include <fparameters\parameters.h>
#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>
#include <fmath\physics.h>

#include <boost/property_tree/ptree.hpp>

int jetAGN()
{
	std::string folder{ prepareOutputfolder() };

	try {
		GlobalConfig = readConfig();
		prepareGlobalCfg();
		State model(GlobalConfig.get_child("model"));

	//	model.photon.injection.ps.iterate([&model](const SpaceIterator& i){
		//	double nph = blackBody(i.val(DIM_E), i.val(DIM_R));
			//model.nph.set(i, nph);
		//});


	//	radiativeLosses(model);


		ParamSpaceValues psv(model.electron.ps);

		psv.fill([&](const SpaceIterator& i){
			double E = i.val(DIM_E);
			double z = i.val(DIM_R);
			return frad(E, z);
		});

		//writeEmax(folder + "\\Emax", model.electron);
		//writeEandRParamSpace(folder + "\\frad", psv, model.electron.ps[DIM_T].size()-1);

		injection(model.electron, model);
		

		std::cout << "checking injected power:" << '\t' << computeInjectedPower(model.electron.injection, 0) << std::endl;		

	//	writeAllSpaceParam("electronInj.txt", model.electron.injection);
		writeEnergyFunction(folder+"\\electronInj_E_r0", model.electron.injection, 0, 0); //escribe Q(E), para r(0) y t(0)
		writeEnergyFunction(folder+"\\electronInj_E_rmax", model.electron.injection, model.electron.ps[1].size()-1, 0); //escribe Q(E), para r(0) y t(0)
	
		distribution(model.electron, model);
	
		writeAllSpaceParam(getFileName(folder, "electronDist"), model.electron.distribution);
		writeEandTParamSpace(getFileName(folder, "electronDist_ET"), model.electron.distribution, model.electron.ps[1].size() / 2);
		writeRandTParamSpace(getFileName(folder,"electronDist_RT"), model.electron.distribution, model.electron.ps[0].size()/2);
		writeEnt(getFileName(folder, "E_NT_r"), model.electron.distribution);


		//lo siguiente es una funcion rapida para llenar N(E) asi pruebo las luminosidades
		//model.electron.distribution.fill([&model](const SpaceIterator& i){
		//	return model.electron.injection.get(i);
		//});

		//processes(model, getFileName(cfg,folder,"luminosity"));

	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
		throw;
	}

	return 0;
}

