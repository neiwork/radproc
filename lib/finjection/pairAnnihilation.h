#pragma once

#include <fparticle\Particle.h>

class AnnihilationData {
public:
	double  E;
	Vector Nep; //esto es para los electrones primarios
	Vector Eep;
	Vector Ne;
	Vector Ee;
	Vector Np; //la p es de positron
	Vector Ep;

};

/* pairAnnihilation calculates the photon injection due to pair annihilation -> Coppi & Blandford 1990 esta esta mal! */ 
double pairAnnihilation(double E, Particle& electron, Particle& secondaryElectron, Particle& positron);

/* pairAnnihilation calculates the photon injection due to pair annihilation -> Boettcher & Schlickeiser 1996
                                                                                     Svensson 1982*/ 