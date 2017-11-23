#pragma once

#include <fparticle\particle.h>

/*Atoyan & Dermer 2003; la seccion eficaz e inelasticidad del canal de pares es de Begelman, Rudak & Sikora 1990*/ 
//double lossesPhotoHadronic(double E, Particle& particle, const ParamSpaceValues& tpf, 
//	                     const SpaceCoord& psc, double phEmin, double phEmax);

double lossesPhotoHadronic(double E, Particle& particle, fun1 tpf,
							double phEmin, double phEmax);


//estos los comparto porque los necesito para la luminosidad
double cPionPH(double u);   //limite inferior
double cPairPH(double u);   //limite inferior

double dPH(double u, double E, double mass); //limite superior