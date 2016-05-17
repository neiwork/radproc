#pragma once

#include <fmath/interpolation.h>

class ParamSpaceValues;
class Dimension {
public:
	
	Vector values; // all the values in the discretization of this dimension

	inline double first() const {
		return values.front();
	};

	inline double last() const {
		return values.back();
	};

	Dimension(size_t size, std::function<void(Vector&)> initializer);

	//void update(Parameters& pars, size_t index);

	const double& operator[](const size_t& index) const;

	size_t size() const;
};

