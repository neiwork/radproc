#pragma once

#include <fmath\physics.h>
#include <fparameters\ParamSpace.h>

enum ParticleType { PT_proton, 
                    PT_electron, 
                    PT_photon, 
                    PT_positron, 
                    PT_pion, 
                    PT_muon,
					PT_muon_R,  //incluye el muR- y el muL+
					PT_muon_L,  //incluye el muL- y el muR+
					PT_muon_L_plus,
					PT_muon_L_minus,
					PT_muon_R_plus,
					PT_muon_R_minus,
					PT_neutrino,
					PT_muonNeutrino,
					PT_muonAntiNeutrino,
					PT_electronNeutrino,
					PT_electronAntiNeutrino,
                    PT_pair,
					PT_secondaryElectron,
					PT_neutron,
                    PT_dummy };

/*
	Particle State
	
	Among other things, it has the current distribution of particles of a given type within the model.
*/
class Particle {
public:

	ParticleType type;

	double mass;

	double logEmax;
	double logEmin;

	double emax() const;
	double emin() const;

	static void initializeEnergyPoints( Vector& energyPoints, double logEmin, double logEmax );

	Particle(ParticleType t, double m, double emin, double emax, int nE);

	
	/* Creates the vectors for injection and distribution 
	   according to the registered dimensions. */
	void initialize();

	ParamSpace ps;

	ParamSpaceValues injection;
	ParamSpaceValues distribution;

	Dimension* eDim() const;
	
	fun1 eInterpolator(ParamSpaceValues psv) const;
	
	//void iterate(std::function <void(const SpaceIterator&)> b); 
	//void iterate(std::function <void(const SpaceIterator&)> b, std::initializer_list<int> fixedDimensions = {}) const;

	/* Interpolates and returns the distribution of a particle at a given energy level */
	double dist(double e) const;


};

class Proton : public Particle {
public:
	Proton() :Particle(PT_proton, protonMass, protonLogEmin, protonLogEmax, nEnergies) {

	}
};

class Electron : public Particle {
public:
	Electron() : Particle(PT_electron, electronMass, electronLogEmin, electronLogEmax, nEnergies) {
	
	}
};

class Photon : public Particle {
public:
	Photon() : Particle(PT_photon, 0.0, photonLogEmin, photonLogEmax, nPhotonEnergies) {
	}
};

class Pion : public Particle {
public:
	Pion() : Particle(PT_pion, chargedPionMass, pionLogEmin, pionLogEmax, nEnergies){}

};

class Muon : public Particle {
public:
	Muon() : Particle(PT_muon, muonMass, muonLogEmin, muonLogEmax, nEnergies){}

};

class Neutrino : public Particle {

};

class Neutron : public Particle {

};

class SecondaryElectron : public Particle {
public:
	SecondaryElectron() :Particle(PT_secondaryElectron, electronMass, electronLogEmin, 12.0, nEnergies) {}
};

class Positron: public Particle {
public:
	Positron() : Particle(PT_positron, electronMass, electronLogEmin, 12.0, nEnergies) {}   //estas particulas son los positrones
};

