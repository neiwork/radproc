#include <stdio.h>


#include "maximumEnergy.h"
#include "checkPower.h"
#include "processes.h"
#include "radiativeLosses.h"
#include "nonThermalLuminosity.h"
#include "distribution.h"
#include "injection.h"
#include "oneZone.h"
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
			
		//radiativeLosses(model,folder+"\\electronLosses.txt");
		
		//frad(z)
		static const double Gamma = GlobalConfig.get<double>("Gamma");
		static const double starT = GlobalConfig.get<double>("starT");

		double E = P2(electronMass*cLight2) / (boltzmann*starT) / Gamma;
		std::ofstream file;
		
		file.open("frad.txt", std::ios::out);

		model.electron.ps.iterate([&](const SpaceIterator& i) {
			double z = i.val(DIM_R);
			double f = frad(E, z);
			file << z/pc << '\t' << f << std::endl;
		}, { 0, -1, 0 });

		//ParamSpaceValues psv(model.electron.ps);

		//psv.fill([&](const SpaceIterator& i){
		//	double E = i.val(DIM_E);
		//	double z = i.val(DIM_R);
		//	return frad(E, z);
		//});

		//writeEandRParamSpace(folder + "\\frad", psv, model.electron.ps[DIM_T].size()-1);
		//writeEmax(folder + "\\Emax", model.electron);

		injection(model.electron, model);
		
		double totalL = computeInjectedPower(model.electron.injection, 0);

		std::cout << "checking injected power:" << '\t' << totalL << std::endl;		

		//writeAllSpaceParam(folder+"\\electronInj.txt", model.electron.injection);
		writeEnergyFunction(folder+"\\electronInj_E_r0", model.electron.injection, 0, 0); //escribe Q(E), para r(0) y t(0)
		writeEnergyFunction(folder+"\\electronInj_E_rmax", model.electron.injection, model.electron.ps[1].size()-1, 0); //escribe Q(E), para rmax y t(0)
	
		//oneZoneDistribution(model.electron, model);
		//processesOneZone(model, folder + "\\oneZoneLum.txt");

		distribution(model.electron, model);
		//distWOLosses(model.electron, model);


		//static const double Gamma = GlobalConfig.get<double>("Gamma");
		double totalE = computeInjectedEnergy(model.electron.distribution); 
		double totalE2 = totalL*(model.electron.ps[1].last() - model.electron.ps[1].first()) / (cLight*Gamma);
		std::cout << "checking total energy:" << '\t' << totalE << '\t' << totalE2 << std::endl;

	
		//writeAllSpaceParam(getFileName(folder, "eDist"), model.electron.distribution);
		//writeEandTParamSpace(getFileName(folder, "eDist_ET"), model.electron.distribution, model.electron.ps[1].size() / 2);
		//writeRandTParamSpace(getFileName(folder, "eDist_RT"), model.electron.distribution, model.electron.ps[0].size() / 2);
		writeEandRParamSpace(getFileName(folder, "eDist_final"), model.electron.distribution, model.electron.ps[2].size() - 1);
		writeEnt(getFileName(folder, "E_NT_r"), model.electron.distribution);


		processes(model, getFileName(folder, "luminosity"));

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