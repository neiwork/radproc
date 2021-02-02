#pragma once

#include <fmath/physics.h>
#include <fparticle/Particle.h>

/* Reynoso, Medina & Romero 2011, A&A; Rybicki & Lightman 1979
Factor to account for the effect of synchrotron self absorption (SSA) within the corona */ 

double fSSA(double x, double E, const Particle& p, const double magneticField, const SpaceCoord& psc);
double fSSA2(double x, double E, const Particle& p, const double magneticField, const SpaceCoord& psc);

/*esta rutina esta implementada en el proyecto porque depende de las coordenadas*/
//double opticalDepthSSA(double E, double mass, double Emin, double Emax, const Particle& creator);
