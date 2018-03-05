#pragma once

#include <fparticle\Particle.h>



double dLnt(double z, double Gc, double z_int, double Rs);

double frad(double E, double z, double Gc);// , double Grel);

double frad_2(double E, double z, double Gc, Particle& electron);
