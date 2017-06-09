#include "jetSN.h"


#include <stdio.h>



#include "checkPower.h"
#include "processes.h"
#include "dynamics.h"
#include "radiativeLosses.h"
#include "nonThermalLuminosity.h"
#include "distribution.h"
#include "injection.h"
//#include "oneZone.h"
#include "state.h"
#include "modelParameters.h"

#include "write.h"
//#include "checks.h"

#include <inout\ioutil.h>
#include <fparticle\Particle.h>
#include <fparameters\parameters.h>
#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>
#include <fmath\physics.h>

#include <boost/property_tree/ptree.hpp>

int jetSN()
{
	std::string folder{ prepareOutputfolder() };

	try {
		GlobalConfig = readConfig();
		prepareGlobalCfg();
		State model(GlobalConfig.get_child("model"));
			
		//radiativeLosses(model,folder+"\\electronLosses.txt", Gc);
				
		int nR = model.electron.ps[DIM_R].size();// -1;

		Vector Gc(nR, 0.0); //ver como inicializarlo
		
		gammaC(model, Gc);

		writeEvol(folder + "\\Gc", model.electron.ps, Gc);

		injection(model.electron, model, Gc);
		writeAllSpaceParam(folder + "\\eInj", model.electron.injection);
		
		//double totalL = computeInjectedPower(model.electron.injection, 0);

		//std::cout << "checking injected power:" << '\t' << totalL << std::endl;		

		distribution(model.electron, model, Gc);
		writeAllSpaceParam(folder + "\\eDist", model.electron.distribution);
		
		//static const double Gamma = GlobalConfig.get<double>("Gamma");
		//double totalE = computeInjectedEnergy(model.electron.distribution); 
		
		//std::cout << "checking total energy:" << '\t' << totalE << std::endl;

		//writeAllSpaceParam(getFileName(folder, "eDist"), model.electron.distribution);

		processes(model, folder + "\\luminosity", Gc);

	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
		throw;
	}

	return 0;
}








//lo siguiente es una funcion rapida para llenar N(E) asi pruebo las luminosidades
//model.electron.distribution.fill([&model](const SpaceIterator& i){
//	return model.electron.injection.get(i);
//});