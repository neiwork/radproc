#pragma once

#include "modelParameters.h"
#include <fparticle/Particle.h>
#include <fmath\mathematics.h>

class State {
public:
	std::vector<Particle*> particles;

	Electron electron;
	//Proton proton;
	//Pion pion;
	//Muon muon;
	Photon photon;
	//SecondaryElectron secondaryElectron;
	//Positron positron;

	ParamSpaceValues nph;
	//fun1 tpf;

	State();
	
	void initializeParticle(Particle& p);
};
