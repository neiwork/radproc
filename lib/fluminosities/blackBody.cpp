#include "blackBody.h"
#include <fmath/physics.h>

double bb_RJ(double frequency, double temp) {
	return 2.0*frequency*frequency*boltzmann*temp / cLight2;
}

double bb(double frequency, double temp) {	
	return 2.0*planck*frequency*frequency*frequency / cLight2 / (exp(planck*frequency/(boltzmann*temp))-1.0);
}