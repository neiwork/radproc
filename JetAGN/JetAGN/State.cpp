#include "State.h"

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
	p.ps.add(new Dimension(nR + 1, &Parameters::R, bind(initializeRPoints, _1, rmin, rmax)));
	p.ps.add(new Dimension(nTimes + 1, &Parameters::T, bind(initializeLinearPoints, _1, timeMin, timeMax)));

	p.ps.addDerivation([](const SpaceIterator& i){
		derive_parameters_r(i.par.E, i.par.R, i.par.T);
	});

	p.initialize();
}




