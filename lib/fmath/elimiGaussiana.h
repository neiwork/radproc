#pragma once

#include <iostream>
#include "physics.h"

/** elimiGaussiana solves a system of n equation,
   with an extended matrix coefficient a(n x n+1).
	 @warning a is modified in the process
	 @returns result in Ne*/ 
void elimiGaussiana (int N2, Matrix& a, Vector& Ne);

