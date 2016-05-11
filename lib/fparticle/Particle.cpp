#include "Particle.h"

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

	ps.add(new Dimension(nE + 1, &Parameters::E, bind(initializeEnergyPoints, _1, emin, emax)));
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

double Particle::dist(double e) const
{
	return eDim()->interpolate(e, distribution);
}

//(std::function<void(const SpaceIterator&)> body, std::initializer_list<int> fixedDimensions) const
//{
	//iterate(Parameters(parameters), body, fixedDimensions);
//void Particle::iterate(std::function <void(const SpaceIterator&)> b, std::initializer_list<int> fixedDimensions) const
//{
//	ps.iterate(b, fixedDimensions);
//}

Dimension* Particle::eDim() const
{
	return ps.dimensions[0];
}

fun1 Particle::eInterpolator(ParamSpaceValues psv) const
{
	return psv.dimInterpolator(0);
}

double Particle::emin() const
{
	return 1.6e-12*pow(10.0, logEmin);
}

double Particle::emax() const
{
	return 1.6e-12*pow(10.0, logEmax);
}

SecondaryElectron::SecondaryElectron() :Particle(PT_secondaryElectron, electronMass, parameters.electronLogEmin, 12.0, parameters.nEnergies)
{

}

Positron::Positron() : Particle(PT_positron, electronMass, parameters.electronLogEmin, 12.0, parameters.nEnergies)
{

}

Muon::Muon() : Particle(PT_muon, muonMass, parameters.muonLogEmin, parameters.muonLogEmax, parameters.nEnergies)
{

}

Proton::Proton() : Particle(PT_proton, protonMass, parameters.protonLogEmin, parameters.protonLogEmax, parameters.nEnergies)
{

}

Electron::Electron() : Particle(PT_electron, electronMass, parameters.electronLogEmin, parameters.electronLogEmax, parameters.nEnergies)
{

}

Photon::Photon() : Particle(PT_photon, 0.0, parameters.photonLogEmin, parameters.photonLogEmax, parameters.nPhotonEnergies)
{

}

Pion::Pion() : Particle(PT_pion, chargedPionMass, parameters.pionLogEmin, parameters.pionLogEmax, parameters.nEnergies)
{

}
