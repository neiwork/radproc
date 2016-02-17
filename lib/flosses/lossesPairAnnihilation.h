#pragma once


#include <fparticle\Particle.h>


/*This function calculates the cooling time for positrons due to pair annihilation*/
//double lossesPairAnnihilation(double E, void* data, Particle& particle);
double lossesPairAnnihilation(double E, Particle& particle, Particle& secondaryParticle);

/*This function calculates the cooling time for electrons due to pair annihilation*/
//double LossesAnnEl(double E, const CoupledEqSys* previous, CoupledEqSys* current);

