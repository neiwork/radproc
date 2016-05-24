#pragma once

#include <fparticle\particle.h>

#include <string>
#include <iostream>
#include <fstream>


/* lossesAnisotropicIC calculates the energy loss [erg/s] due to anisotropic IC*/
double lossesAnisotropicIC(double E, Particle& particle, double r);

/*these three functions are used by luminosityAniIC*/
double b_theta(double theta, double w0, double E);

//double difN(double theta, double w, double w0, double E, double r);

class IntLossesOpt {
public:
	int samples_x;
	int samples_t;
	int samples_y;
};

extern IntLossesOpt DefOpt_IntLosses;
double intTriple(double E, double eps_min, double eps_max, double r, fun2 c, fun2 d, fun3 f, const IntLossesOpt& opt = DefOpt_IntLosses);
