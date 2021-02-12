#pragma once

#include <fmath\mathematics.h>
#include <fparticle\Particle.h>
#include <boost/property_tree/ptree_fwd.hpp>

//
//class Electron : public ParticleCfg<Electron> {};
//class Photon : public ParticleCfg<Photon> {};

const DimensionCoord
	DIM_E = 0;

/* define the inital values of the global parameters*/
void prepareGlobalCfg();

void initializeEnergyPoints(Vector& v, double logEmin, double logEmax);

double computeMagField(double z);

double computeDens(double z);

double jetRadius(double z, double openingAngle);

double eEmax(double mass, double z, double B);


