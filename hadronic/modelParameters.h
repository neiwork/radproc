#pragma once

#include <fmath\mathematics.h>
#include <fparticle\Particle.h>
#include <boost\property_tree\ptree_fwd.hpp>

//
//class Electron : public ParticleCfg<Electron> {};
//class Photon : public ParticleCfg<Photon> {};

const DimensionCoord
	DIM_E = 0,
	DIM_R = 1,
	DIM_T = 2;

/* define the inital values of the global parameters*/
void prepareGlobalCfg();

void initializeRPoints(Vector& v, double Rmin, double Rmax);

void initializeCrossingTimePoints(Vector& time, double rMin, double rMax);

void initializeEnergyPoints(Vector& v, double logEmin, double logEmax);

/* change_parameters changes the values of some parameters for the iteration on variable r*/
//void derive_parameters_r(double E, double R, double T);

//double stagnationPoint(double z);

double beta(double gamma);

double computeMagField(double z);

double zRecol();

double jetRadius(double z);

double eEmax(double z, double B);

double stagnationPoint(double z);


//double vWind(double r, double starR);

