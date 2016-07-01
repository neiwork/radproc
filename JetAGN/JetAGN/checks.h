#pragma once

#include "State.h"

bool check_vec(const std::vector<double>& a, const std::vector<double>& b);

void check_inj(const Particle& p);

void check_dist(const Particle& p);

void check_refactor(const Particle& p);

void testCanPeek();

void test3DPrint();

void testBinarySearch();

void testMultiDimensionalInterpolate();

void zeroToN(Vector& v);


void kolmogorov();