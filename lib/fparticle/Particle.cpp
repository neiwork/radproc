#include "Particle.h"

#include <fparameters/Dimension.h>

namespace {
	void zeroToN(Vector& v) {
		for (int i = 0; i < (int)v.size(); ++i) {
			v[i] = i;
		}
	}
}
//si la inicializo asi entonces el modelo NO depende del tiempo
Particle::Particle(ParticleType t, double m, double emin, double emax, int nE)
	:injection(ps,false), // these PSVs are not initialized immediately, only after this PS has been constructed
	distribution(ps,false)
{
	using std::bind; using namespace std::placeholders; // para _1, _2, etc.

	type = t;
	mass = m;
	logEmin = emin;
	logEmax = emax;

	ps.add(new Dimension(nE + 1, bind(initializeEnergyPoints, _1, emin, emax)));
}
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

void Particle::initialize() {
	injection.initialize();
	distribution.initialize();
}

void Particle::initializeEnergyPoints(Vector& v, double logEmin, double logEmax)
{
	double Emax  = 1.6e-12*pow(10,logEmax);    
	double Emin  = 1.6e-12*pow(10,logEmin);

	double E_int = pow((10*Emax/Emin),(1.0/(v.size()-1)));

	v[0] = Emin;

	for (size_t i=1; i < v.size() ; ++i){  
		v[i] = v[i-1]*E_int;
	}

}


//
//void Particle::initializeLinearEnergyPoints( int n )
//{
//	double Emax  = 1.6e-12*pow(10,logEmax);    
//	double Emin  = 1.6e-12*pow(10,logEmin);
//
//	double E_int = (Emax-Emin)/n;
//
//	energyPoints[0] = Emin;
//
//	for (size_t i=1; i < energyPoints.size() ; ++i){  
//		energyPoints[i] = energyPoints[i-1]+E_int;
//	}
//
//}

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
ParticleConfig ParticleCfg<Proton>::config{ PT_proton, protonMass, 0, 0, 0 };
ParticleConfig ParticleCfg<Electron>::config{ PT_electron, electronMass, 0, 0, 0 };
ParticleConfig ParticleCfg<Photon>::config{ PT_photon, 0, 0, 0, 0 };
ParticleConfig ParticleCfg<Pion>::config{ PT_pion, 0, 0, 0, 0 };
ParticleConfig ParticleCfg<Muon>::config{ PT_muon, muonMass, 0, 0, 0 };
ParticleConfig ParticleCfg<Neutrino>::config{ PT_neutrino, 0, 0, 0, 0 };
ParticleConfig ParticleCfg<Neutron>::config{ PT_neutron, neutronMass, 0, 0, 0 };
ParticleConfig ParticleCfg<SecondaryElectron>::config{ PT_secondaryElectron, 0, 0, 0, 0 };
ParticleConfig ParticleCfg<Positron>::config{ PT_positron, 0, 0, 0, 0 };