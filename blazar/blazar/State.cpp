#include "State.h"

#include "modelParameters.h"

#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\parameters.h>

#include <boost/property_tree/ptree.hpp>

State::State(boost::property_tree::ptree& cfg) :
 electron{ "electron" },
 proton{ "proton" },
 photon{ "photon" },
 tpf(photon.ps, false),
 magf(photon.ps, false)
 {
	static const double Mbh = GlobalConfig.get<double>("Mbh")*solarMass;
	static const double zInt = GlobalConfig.get<double>("Rdiss");
	 
	 particles.push_back(&electron);
	 particles.push_back(&proton);
	particles.push_back(&photon);

	for (auto p : particles) {
		initializeParticle(*p,cfg);
	}

	magf.initialize(computeMagField(zInt));
}


Dimension* State::createDimension(Particle& p, string dimid,
	function<void(Vector&, double, double)> initializer, function<double(double)> to_linear, function<double(double)> from_linear, boost::property_tree::ptree& cfg)
{
	int samples = p.getpar<int>(cfg, "dim." + dimid + ".samples");
	double min = p.getpar<double>(cfg, "dim." + dimid + ".min");
	double max = p.getpar<double>(cfg, "dim." + dimid + ".max");
	return new Dimension(samples, bind(initializer, placeholders::_1, min, max), to_linear, from_linear);
}

auto l10 = [](double x) { return (x > 0.0) ? log10(x) : -300.0; };
auto e10 = [](double x) { return std::pow(10.0,x); };
/*Dimension* State::createDimension(Particle& p, std::string dimid, std::function<void(Vector&,double,double)> initializer, boost::property_tree::ptree& cfg) {
	int samples = p.getpar<int>(cfg,"dim." + dimid + ".samples");
	double min = p.getpar<double>(cfg, "dim." + dimid + ".min");
	double max = p.getpar<double>(cfg, "dim." + dimid + ".max");
	return new Dimension(samples, bind(initializer, std::placeholders::_1, min, max));
}*/

void State::initializeParticle(Particle& p, boost::property_tree::ptree& cfg)
{
	using std::bind;

	p.configure(cfg.get_child("particle.default"));
	p.configure(cfg.get_child("particle."+p.id));

	p.ps.add(createDimension(p,"energy",initializeEnergyPoints, l10, e10, cfg));

	// we can't use createDimension because we're multiplying by pc before creating them
	// add dimension for R
	/*
	double rmin = p.getpar(cfg,"dim.radius.min", 1.0)*pc;
	double rmax = p.getpar(cfg,"dim.radius.max", 1.0e3)*pc;
	int nR = p.getpar(cfg,"dim.radius.samples", 5); // solo por ahora; y no deberia ser usado directamente desde otro lado
	p.ps.add(new Dimension(nR, bind(initializeRPoints, std::placeholders::_1, rmin, rmax)));

	// add dimension for T
	double tmin = p.getpar(cfg, "dim.time.min", 1.0)*pc;
	double tmax = p.getpar(cfg, "dim.time.max", 1.0e3)*pc;
	int tR = p.getpar(cfg, "dim.time.samples", 5); // solo por ahora; y no deberia ser usado directamente desde otro lado
	p.ps.add(new Dimension(tR, bind(initializeCrossingTimePoints, std::placeholders::_1, tmin, tmax)));
	*/

	//p.ps.addDerivation([](const SpaceIterator& i){
	//	derive_parameters_r(i.val(DIM_E), i.val(DIM_R), i.val(DIM_T));
	//});

	p.initialize();
}




