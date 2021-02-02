#pragma once

#include <fparticle/Particle.h>

/*pion injection due to photoHadronic interaction
(Atoyan & Dermer 2003; Vila & Aharonian 2009) */

double pgammaPionInj(double E, const Particle& creator, const ParamSpaceValues& tpf, const SpaceCoord& psc, 
						double tpEmin, double tpEmax);
//double pgammaPionInj(double E, const Particle& creator,	const ParamSpaceValues& tpf, const SpaceCoord& distCoord);

//estas funciones las muestro porque las use para la luminosidad
//aca solo incluyo la produccion de piones

//double t_pion_PH(double E, const Particle& particle, fun1 tpf, double tpEmin, double tpEmax);   

//double omegaPH(double E, const Particle& particle, fun1 tpf, double tpEmin, double tpEmax);


double omegaPH(double E, const Particle& particle, const ParamSpaceValues& tpf, const SpaceCoord& distCoord,
					double tpEmin, double tpEmax);  //E=Ep  
//double omegaPHsimple(double E, const Particle& particle, const ParamSpaceValues& tpf, const SpaceCoord& distCoord,
//					double tpEmin, double tpEmax);  //E=Ep  

double t_pion_PH(double E, const Particle& particle, const ParamSpaceValues& tpf, const SpaceCoord& distCoord,
					double tpEmin, double tpEmax);
//double t_pion_PHsimple(double E, const Particle& particle, const ParamSpaceValues& tpf, const SpaceCoord& distCoord,
//					double tpEmin, double tpEmax);