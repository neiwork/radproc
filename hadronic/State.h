#pragma once

#include "modelParameters.h"
#include <fparticle/Particle.h>
#include <fmath\mathematics.h>
#include <fmath\physics.h>
#include <boost/property_tree/ptree_fwd.hpp>

class State {
public:
	std::vector<Particle*> particles;

	Particle proton;
	Particle photon;

	//ParamSpaceValues magf;

	State(boost::property_tree::ptree& cfg);
	
	static Dimension* createDimension(Particle& p, std::string dimid, std::function<void(Vector&, double, double)> initializer, boost::property_tree::ptree& cfg);

	static void initializeParticle(Particle& p, boost::property_tree::ptree& cfg);
};
