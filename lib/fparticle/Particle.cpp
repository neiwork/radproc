#include "Particle.h"

#include <fparameters/Dimension.h>
#include <boost/property_tree/ptree.hpp>

/*template<typename T> T Particle::getpar(boost::property_tree::ptree& cfg, const std::string &path, const T& def ) {
	return cfg.get<T>("particle." + id + "." + path, cfg.get<T>("particle.default." + path, def));
}


template int Particle::getpar<int>(boost::property_tree::ptree& cfg, const std::string &path, const int& def );
template double Particle::getpar<double>(boost::property_tree::ptree& cfg, const std::string &path, const double& def );
*/  //version linux

namespace {
	void zeroToN(Vector& v) {
		for (int i = 0; i < (int)v.size(); ++i) {
			v[i] = i;
		}
	}
}


Particle::Particle(const std::string& id)
:id(id),
 mass{ 0 },
injection{ ps, false }, // these PSVs are not initialized immediately, only after this PS has been constructed
 distribution{ ps, false }
{
}


Particle::Particle(Particle&& other):
	ps(other.ps),
	injection(this->ps, false),
	distribution(this->ps, false),
	id(other.id),
	mass(other.mass),
	logEmax(other.logEmax),
	logEmin(other.logEmin)
{
	injection.values = std::move(other.injection.values);
	distribution.values = std::move(other.distribution.values);
}




void Particle::configure(boost::property_tree::ptree& cfg) {
	mass = cfg.get<double>("mass", mass);
	logEmin = cfg.get<double>("dim.energy.min", logEmin);
	logEmax= cfg.get<double>("dim.energy.max", logEmax);
}

void Particle::initialize() {
	injection.initialize();
	distribution.initialize();
}


Dimension* Particle::eDim() const
{
	return ps.dimensions[0];
}

double Particle::emin() const
{
	return 1.6e-12*pow(10.0, logEmin);
}

double Particle::emax() const
{
	return 1.6e-12*pow(10.0, logEmax);
}



////si la inicializo asi entonces el modelo NO depende del tiempo
//Particle::Particle(ParticleType t, double m, double emin, double emax)
//	:injection(ps,false), // these PSVs are not initialized immediately, only after this PS has been constructed
//	distribution(ps,false)
//{
//	type = t;
//	mass = m;
//	logEmin = emin;
//	logEmax = emax;
//}


//
//
////si la inicializo asi entonces el modelo depende del tiempo  (o es inhomogeneo?? )
//Particle::Particle(ParticleType t, double m, double emin, double emax, int nE, int nT, double tMin, double tStep)
//	:energyPoints(nE+1, 0.0),  //el +1 es por el paso logaritmico
//	timePoints(nT+1, 0.0),
//	injection( energyPoints.size()*timePoints.size(), 0.0 ),
//	distribution( energyPoints.size()*timePoints.size(), 0.0 ),
//	injection(ps),
//	distribution(ps)
//{
//	type = t;
//	mass = m;
//	logEmin = emin;
//	logEmax = emax;
//	
//	initializeEnergyPoints(energyPoints,logEmin,logEmax);
//	initializeTimePoints(nT,tMin,tStep); 
//
//	ps.add(new Dimension(nE + 1, &Parameters::E, [emin,emax](Vector& v){
//		initializeEnergyPoints(v, emin, emax);
//	}));
//
//
//	initialize();
//}
