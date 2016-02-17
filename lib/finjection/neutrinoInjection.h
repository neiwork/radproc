#pragma once


#include <fparticle\particle.h>


/* calculates muonNeutrino by decay of pion+ and muon- (L y R)  (Lipari 2007) */ 
double muonNeuInj(double E, int pos_E, int pos_t, Particle& creator);

/* calculates muonAntiNeutrino by decay of pion- and muon+ (L y R) (Lipari 2007) */ 
double muonAntiNeuInj(double E, int pos_E, int pos_t, Particle& creator);

/* calculates electronNeutrino by decay of mu+ (L y R) (Lipari 2007) */ 
double electronNeuInj(double E, int pos_E, int pos_t, Particle& creator);

/* calculates electronAntiNeutrino by decay of mu- (L y R) (Lipari 2007) */ 
double electronAntiNeuInj(double E, int pos_E, int pos_t, Particle& creator);