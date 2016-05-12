#pragma once

#include <fmath/interpolation.h>
#include "parameters.h"

class ParamSpaceValues;
class Dimension {
public:
	
	Vector values; // all the values in the discretization of this dimension
	
	double min;
	double max;
	
	double Parameters::*fieldRef; // the pointer-to-member used to update the parameters object
	
	Dimension(size_t size, std::function<void(Vector&)> initializer, double Parameters::*fieldRef = nullptr);

	void update(Parameters& pars, size_t index);

	double interpolate(double x, const ParamSpaceValues& psv) const;

	const double& operator[](const size_t& index) const;

	size_t size() const;
};

