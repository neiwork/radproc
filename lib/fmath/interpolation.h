#pragma once

#include "physics.h"

//double interpolMod(double& E, const Vector& ener, const Vector& lum, const int last);


//double interpol(double& E, const Vector& ener, const Vector& lum, const int last);
double interpol(double& E, const Vector& ener, const Vector& lum, const int last, const int first=0);
double interpolNew(double& E, const Vector& key, const Vector& val, const int size, const int first = 0);  //esta es con binary search

double interpolDoble(double E, double t, const Vector& ener, const Vector& time, const Vector& dist);
