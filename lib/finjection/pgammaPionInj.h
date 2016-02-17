#pragma once

#include <fparticle\particle.h>

/*pion injection due to photoHadronic interaction
(Atoyan & Dermer 2003; Vila & Aharonian 2009) */
double pgammaPionInj(double E, Vector Nproton, Particle& particle, Particle& proton, fun1 tpf);

//double pgammaPairInj(double E, Vector Nproton, Particle& particle, Particle& proton, fun1 tpf)  ;


//estas funciones las muestro porque las use para la luminosidad
double t_pion_PH(double E, Particle& particle, fun1 tpf);   //aca solo incluyo la produccion de piones

double omegaPH(double E, Particle& particle, fun1 tpf);

//double fOmegaPHPion(double u, double t, double E, double mass, fun1 tpf);
//
//double f_t_PHPion(double u, double t, fun1 tpf)


