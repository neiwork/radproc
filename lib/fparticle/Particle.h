#pragma once

#include <fmath\physics.h>

#include <fparameters\ParamSpaceValues.h>
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

struct ParticleConfig {
public:
	ParticleType type;
	double mass;
	double logEmin;
	double logEmax;
	int nE;
};

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
	Particle(const ParticleConfig& pcfg) :Particle(pcfg.type, pcfg.mass, pcfg.logEmin, pcfg.logEmax, pcfg.nE){};

	
	/* Creates the vectors for injection and distribution 
	   according to the registered dimensions. */
	void initialize();

	ParamSpace ps;

	ParamSpaceValues injection;
	ParamSpaceValues distribution;

	Dimension* eDim() const;
};


template <class ConfigHolder>
class ParticleCfg : public Particle {
public:
	static ParticleConfig config;
	ParticleCfg():Particle(ParticleCfg::config) {}
};