#include "blazar.h"
#include "ioutil.h"

#include "chiSquare.h"
#include "write.h"
#include "photonInjection.h"
#include "photonInjection.h"
#include "distribution.h"
#include "injection.h"
#include "State.h"
#include "modelParameters.h"

#include <fparameters/parameters.h>
#include <fparticle/Particle.h>

#include <boost/property_tree/ptree.hpp>





void blazar() {

	std::string folder{ prepareOutputfolder() };

	Matrix data;
	readData("BlLac.txt", data);

	GlobalConfig = readConfig();
	prepareGlobalCfg();
	
	State model(GlobalConfig.get_child("model"));


	//	model.photon.injection.ps.iterate([&model](const SpaceIterator& i){
	//	double nph = blackBody(i.val(DIM_E), i.val(DIM_R));
	//model.nph.set(i, nph);
	//});


	//	radiativeLosses(model);

	injection(model.electron, model);

	distribution(model.electron, model);

	photonTarget(model.photon, model);

	distribution(model.electron, model);

	photonDistribution(model.photon, model);

	writeAllSpaceParam(getFileName(folder, "electronDist"), model.electron.distribution);
	writeAllSpaceParam(getFileName(folder, "SED"), model.photon.distribution);

	Matrix data;
	readData("PKS0048-097_lum.txt", data);
	double chi = chiSquare(model.photon.distribution, data, 2);

	std::cout << chi;

}