#include "State.h"

#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>

State::State() :
	nph(photon.ps,false)
{
	particles.push_back(&electron);
	//particles.push_back(&proton);
	//particles.push_back(&pion);
	//particles.push_back(&muon);
	particles.push_back(&photon);
	//particles.push_back(&secondaryElectron);
	//particles.push_back(&positron);

	for (auto p : particles) {
		initializeParticle(*p);
	}

	nph.initialize();
	tpf = nph.dimInterpolator(0);
}

void State::initializeParticle(Particle& p)
{
	using std::bind; using namespace std::placeholders; // para _1, _2, etc.

	// add dimension for R
	p.ps.add(new Dimension(parameters.nR + 1, bind(initializeRPoints, _1, parameters.rmin, parameters.rmax), &Parameters::R));
	p.ps.add(new Dimension(parameters.nR + 1, bind(initializeCrossingTimePoints, _1, parameters.rmin, parameters.rmax), &Parameters::T));

	p.ps.addDerivation([](const SpaceIterator& i){
		derive_parameters_r(i.par.E, i.par.R, i.par.T);
	});

	p.initialize();
}




