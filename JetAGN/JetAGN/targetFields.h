#pragma once

#include <fmath\physics.h>

#include <fparticle\particle.h>

double cmbBlackBody(double Ep);

double starDensity(double z);

double starBlackBody(double E, double r);

double starIR(double Ep, double z);

//void targetPhotonEnergies(double& EphminS, double& EphminCMB);

void tpfPSV(ParamSpaceValues& psv, fun2 tpf, Particle photon, double Lorentz);