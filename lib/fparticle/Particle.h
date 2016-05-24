#pragma once

#include <fmath\physics.h>

#include <fparameters\ParamSpaceValues.h>
#include <fparameters\ParamSpace.h>

#include <boost/property_tree/ptree_fwd.hpp>

//struct ParticleConfig {
//public:
//	ParticleType type;
//	double mass;
//	double logEmin;
//	double logEmax;
//	int nE;
//};

/*
	Particle State
	
	Among other things, it has the current distribution of particles of a given type within the model.
*/
class Particle {
public:

	const std::string id;

	double mass;

	double logEmax;
	double logEmin;

	double emax() const;
	double emin() const;

	template<typename T> T getpar(boost::property_tree::ptree& cfg, const std::string &path, const T& def = T{}) {
		return cfg.get<T>("particle." + id + "." + path, cfg.get<T>("particle.default." + path, def));
	}

	Particle(const std::string& id);
	
	/* Creates the vectors for injection and distribution 
	   according to the registered dimensions. */
	void initialize();
	void configure(boost::property_tree::ptree& cfg);

	ParamSpace ps;

	ParamSpaceValues injection;
	ParamSpaceValues distribution;

	Dimension* eDim() const;
};