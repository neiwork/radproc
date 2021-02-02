#include "blazar.h"
#include "ioutil.h"

#include "radiativeLosses.h"
#include "chiSquare.h"
#include "write.h"
#include "photonInjection.h"
#include "photonInjection.h"
#include "distribution.h"
#include "injection.h"
#include "State.h"
#include "modelParameters.h"
#include "tpf.h"

#include <fparameters/parameters.h>
#include <fparticle/Particle.h>

#include <boost/property_tree/ptree.hpp>





void blazar() {

	std::string folder{ prepareOutputfolder() };


	GlobalConfig = readConfig();
	prepareGlobalCfg();
	
	State model(GlobalConfig.get_child("model"));
	

	injection(model.electron, model);
	injection(model.proton, model);

	writeAllSpaceParam(getFileName(folder, "electronInj"), model.electron.injection);

	distribution(model.electron, model);
	
	//ParamSpaceValues Qsyn(model.photon.ps);
	synL(model, model.tpf);

	radiativeLosses(model);

	//distribution(model.electron, model);
	distribution(model.proton, model);

	writeAllSpaceParam(getFileName(folder, "electronDist"), model.electron.distribution);
	writeAllSpaceParam(getFileName(folder, "protonDist"), model.proton.distribution);
	
	//photonDistribution(model.photon, model);

	processes(model, getFileName(folder, "SED"));

	
	
	//dopplerBoost(getFileName(folder, "SED"), model.photon.distribution);
	//writeAllSpaceParam(getFileName(folder, "SED"), model.photon.distribution);


//	static const int height = GlobalConfig.get<int>("height");
//	static const int width = GlobalConfig.get<int>("width");
	
//	Matrix obs;
//	matrixInit(obs, height, width, 0.0);
//	readData("PKS0048-097_lum.txt", obs);
//	double dof = height;
//	double chi = chiSquareFit(model.photon.distribution, obs, dof);

//	std::cout << chi;

}