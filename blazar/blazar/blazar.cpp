#include "blazar.h"
#include "ioutil.h"

#include "modelParameters.h"
#include "State.h"
#include "injection.h"
#include "distribution.h"
#include "photonInjection.h"

#include <fparameters/parameters.h>
#include <fparticle/Particle.h>

#include <boost/property_tree/ptree.hpp>





void blazar() {

	std::string folder{ prepareOutputfolder() };

	GlobalConfig = readConfig();

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

	//escribir

	



}