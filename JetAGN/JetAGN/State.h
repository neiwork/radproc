#pragma once

#include "modelParameters.h"
#include <fparticle/Particle.h>
#include <fmath\mathematics.h>
#include <fmath\physics.h>

class State {
public:
	std::vector<Particle*> particles;

	Particle electron;
	//Proton proton;
	//Pion pion;
	//Muon muon;
	Particle photon;
	//SecondaryElectron secondaryElectron;
	//Positron positron;

	ParamSpaceValues nph;
	//fun1 tpf;

	State(boost::property_tree::ptree& cfg);
	
	void initializeParticle(Particle& p, boost::property_tree::ptree& cfg);
};
