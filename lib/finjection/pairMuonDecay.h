#pragma once

#include <fparticle/Particle.h>


/* muonDecay calculates the pair injection due to muon decay 
 * from Ramaty (1974) */ 
 
//double pairMuonDecay(double E, Particle& particle, Particle& creator) ;
double pairMuonDecay(double E, Particle& c, const SpaceCoord& distCoord);
double pairMuonDecayNew(double E, Particle& c, const SpaceCoord& distCoord);
double pairMuonDecayNew2(double E, Particle& c, const SpaceCoord& distCoord);