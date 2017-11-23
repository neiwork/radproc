#include "State.h"

#include "functions.h"
#include "modelParameters.h"
#include <fparameters/Dimension.h>
#include <fmath/physics.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>

State::State(boost::property_tree::ptree& cfg) :
 electron{ "electron" },
 photon{ "photon" },
 magf(photon.ps, false),
denf(photon.ps,false)
{
	particles.push_back(&electron);
	particles.push_back(&photon);

	for (auto p : particles) {
		initializeParticle(*p, cfg);
	}
	
	magf.initialize();
	magf.fill([&](const SpaceIterator& i){
		static const double beta = GlobalConfig.get<double>("beta");
	    double r = i.val(DIM_R);
	    double theta = i.val(DIM_THETA);
		return sqrt( beta * 24.0 * pi * pressureTot(r, theta) );
	  });
	  
//	tpf.initialize();
//	tpf.fill([&](const SpaceIterator& i){
//		return blackBody(i.val(DIM_E), i.val(DIM_R));
//	  };
	
}

Dimension* State::createDimension(Particle& p, std::string dimid, std::function<void(Vector&,double,double)> initializer, boost::property_tree::ptree& cfg) {
	int samples = p.getpar<int>(cfg,"dim." + dimid + ".samples");
	double min = p.getpar<double>(cfg, "dim." + dimid + ".min");
	double max = p.getpar<double>(cfg, "dim." + dimid + ".max");
	return new Dimension(samples, bind(initializer, std::placeholders::_1, min, max));
}

void State::initializeParticle(Particle& p, boost::property_tree::ptree& cfg)
{
	using std::bind;

	p.configure(cfg.get_child("particle.default"));
	p.configure(cfg.get_child("particle."+p.id));

	p.ps.add(createDimension(p, "energy", initializeEnergyPoints, cfg));

	// we can't use createDimension because we're multiplying by pc before creating them
	// add dimension for R
    double rmin = GlobalConfig.get<double>("rCusp");
    double rmax = 3.0 * GlobalConfig.get<double>("rCenter");
	int nR = p.getpar(cfg,"dim.radius.samples", 5); // solo por ahora; y no deberia ser usado directamente desde otro lado
	p.ps.add(new Dimension(nR, bind(initializePoints, std::placeholders::_1, rmin, rmax)));
    
    // add dimension for theta
    double thetamin =0.0;                           // los defino aca porque no se si puedo poner pi en el .json
    double thetamax = pi/2.0 - 1.e-2;
    int thetaR = p.getpar(cfg, "dim.theta.samples", 5);
    p.ps.add(new Dimension(thetaR, bind(initializePoints, std::placeholders::_1, thetamin, thetamax)));

	// add dimension for T
	// double tmin = p.getpar(cfg, "dim.time.min", 1.0)*pc;
	// double tmax = p.getpar(cfg, "dim.time.max", 1.0e3)*pc;
	// int tR = p.getpar(cfg, "dim.time.samples", 5); // solo por ahora; y no deberia ser usado directamente desde otro lado
	// p.ps.add(new Dimension(tR, bind(initializeCrossingTimePoints, std::placeholders::_1, tmin, tmax)));


	//p.ps.addDerivation([](const SpaceIterator& i){
	//	derive_parameters_r(i.val(DIM_E), i.val(DIM_R), i.val(DIM_T));
	//});

	p.initialize();
}




