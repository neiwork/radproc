#pragma once



#include "State.h"

//double emiToLumi(const ParamSpace& pps, ParamSpaceValues& psv, double E, int t_ix);

void processes(State& st, const std::string& filename, Vector& Gc, Vector& tobs);

double Llab(double Lint);