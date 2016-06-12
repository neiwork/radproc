#include "blazar.h"
#include "ioutil.h"

#include <fparameters/parameters.h>
#include <fparticle/Particle.h>

#include <boost/property_tree/ptree.hpp>

void blazar() {

	std::string folder{ prepareOutputfolder() };

	GlobalConfig = readConfig();

	Particle miParticula("megaPion");

	



}