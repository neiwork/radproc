#pragma once

#include "State.h"



//void blobRadius(State& st, Vector& Gc, Vector& Rc);

double eEmax(double z0, double z, double Gc, double B, double Reff);

//void gammaC(State& st, Vector& Gc);
void gammaC(State& st, Vector& Gc, Vector& Rc, Vector& tobs);

double Fe(double g, double y);

void fillMagnetic(State& st, Vector& Gc);