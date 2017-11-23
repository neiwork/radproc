#pragma once

#include "State.h"



double eEmax(double z, double Gc, double B, double Reff);

//void gammaC(State& st, Vector& Gc);
void gammaC(State& st, Vector& Gc, Vector& Rc, Vector& tobs);

double Fe(double g, double y);

void fillMagnetic(State& st, Vector& Gc);