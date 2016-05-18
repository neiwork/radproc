#include "State.h"

#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>

ParticleConfig ParticleCfg<Electron>::config{ PT_electron, electronMass, 0, 0, 0 };

ParticleConfig ParticleCfg<Photon>::config{ PT_photon, 0, 0, 0, 0 };

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
	//tpf = nph.dimInterpolator(0);
}

void State::initializeParticle(Particle& p)
{
	using std::bind; using namespace std::placeholders; // para _1, _2, etc.

	double rmin = 1.0*pc;
	double rmax = 1.0e3*pc;
	int nR = 4; // solo por ahora; y no deberia ser usado directamente desde otro lado

	// add dimension for R
	p.ps.add(new Dimension(nR + 1, bind(initializeRPoints, _1, rmin, rmax)));
	p.ps.add(new Dimension(nR + 1, bind(initializeCrossingTimePoints, _1, rmin, rmax)));

	p.ps.addDerivation([](const SpaceIterator& i){
		derive_parameters_r(i.val(DIM_E), i.val(DIM_R), i.val(DIM_T));
	});

	p.initialize();
}




