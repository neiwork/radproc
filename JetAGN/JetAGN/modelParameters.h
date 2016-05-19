#pragma once

#include <fmath\mathematics.h>
#include <fparticle\Particle.h>
#include <boost/property_tree/ptree.hpp>

namespace bpt = boost::property_tree;
//
//class Electron : public ParticleCfg<Electron> {};
//class Photon : public ParticleCfg<Photon> {};

const DimensionCoord
	DIM_E = 0,
	DIM_R = 1,
	DIM_T = 2;

/* define the inital values of the global parameters*/
void setParameters( bpt::ptree& parcfg );

void initializeRPoints(Vector& v, double Rmin, double Rmax);

void initializeCrossingTimePoints(Vector& time, double rMin, double rMax);

//void initializeLinearPoints(Vector& v, double tMin, double tMax);

/* change_parameters changes the values of some parameters for the iteration on variable r*/
void derive_parameters_r(double E, double R, double T);

double stagnationPoint(double z);

double jetRadius(double z, double openingAngle);

double eEmax(double z, double B);


//double vWind(double r, double starR);

