#pragma once

#include <fmath\interpolation.h>

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

	Dimension(int size, std::function<void(Vector&)> initializer);

	//void update(Parameters& pars, size_t index);

	const double& operator[](const int& index) const;

	int size() const;
};

