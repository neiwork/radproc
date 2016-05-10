#pragma once

#include <fmath\physics.h>

/*This function integrates the function fun over the volumen of the jet
\int _zmin ^zmax 2*pi*Rj(z)^2 * fun dz  */

double intCilindric(double zMin, double zMax, fun1 fun);
