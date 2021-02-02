#pragma once

#include <fparticle/Particle.h>

//double protonNeutron(double E, Vector Ncreator, Particle& particle, Particle& neutron);

/* function to estimate the injection of p (proton or electron) by neutron decay*/
double injNeutronDecay(double E, Particle& p, Particle& n, const SpaceCoord& distCoord);
double injElectronNeutronDecay(double E, Particle& n, const SpaceCoord& distCoord);