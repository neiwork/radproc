#pragma once

//#include <fmath/interpolation.h>
#include <fmath/physics.h>

class ParamSpaceValues;
class Dimension {
public:
	
	Vector values; // all the values in the discretization of this dimension

	std::function<double(double)> to_linear = [](double x) { return (x > 0.0) ? log10(x) : -300.0; };
	std::function<double(double)> from_linear = [](double x) { return x; };

	inline double first() const {
		return values.front();
	};

	inline double last() const {
		return values.back();
	};
	
	/** Finds an hypotetical index in the values vector where a given dimension value would be,
	 * by interpolation based on the values of a lower bound and upper bound indexes.
	 * 
	 * We can take advantage of knowing the order of the scale and interpolate in a more linear
	 * space.
	 * */
	double interpolateIndex(double value, int lower, int upper) const;
		
	Dimension(int size, std::function<void(Vector&)> initializer);
	Dimension(int size, std::function<void(Vector&)> initializer, std::function<double(double)> to_linear, std::function<double(double)> from_linear);

	//void update(Parameters& pars, size_t index);

	const double& operator[](const int& index) const;

	int size() const;
};

