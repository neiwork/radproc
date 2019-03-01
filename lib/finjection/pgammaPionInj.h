#pragma once

#include <fparticle\particle.h>

/*pion injection due to photoHadronic interaction
(Atoyan & Dermer 2003; Vila & Aharonian 2009) */
double pgammaPionInj(double E, const Particle& creator, const SpaceCoord& psc, fun1 tpf, double tpEmin, double tpEmax);

//double pgammaPairInj(double E, Vector Nproton, Particle& particle, Particle& proton, fun1 tpf)  ;


//estas funciones las muestro porque las use para la luminosidad

//aca solo incluyo la produccion de piones
double t_pion_PH(double E, const Particle& particle, fun1 tpf, double tpEmin, double tpEmax);   

double omegaPH(double E, const Particle& particle, fun1 tpf, double tpEmin, double tpEmax);




