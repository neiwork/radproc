#pragma once

#include "physics.h"

/** gaussianScaledPivoting solves a system of n equation,
with an extended matrix coefficient a(n x n+1).
Partial pivoting: the algorithm selects the entry with largest absolute value 
from the column of the matrix that is currently being considered as the pivot element 
(it is generally sufficient to adequately reduce round-off error).
Scaled partial pivoting: the algorithm selects as the pivot element the entry that is largest 
relative to the entries in its row. This strategy is desirable when entries' large 
differences in magnitude lead to the propagation of round-off error. 
Scaled pivoting should be used in a system like the one below where a row's 
entries vary greatly in magnitude.
@warning a is modified in the process
@returns result in X*/



void gaussianScaledPivoting(int N, Matrix& A, Vector& X);