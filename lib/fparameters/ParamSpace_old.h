#pragma once
#include <vector>
#include <functional>
#include <numeric>
#include <algorithm>
#include <fmath/mathematics.h>
#include "parameters.h"

class ParamSpaceValues;
class Dimension {
public:
	
	Vector values; // all the values in the discretization of this dimension
	
	double min;
	double max;
	
	double Parameters::*par; // the pointer-to-member used to update the parameters object
	
	Dimension(size_t size, double Parameters::*par, std::function<void(Vector&)> initializer);

	void update(Parameters& pars, size_t index);

	double interpolate(double x, const ParamSpaceValues& psv) const;

	size_t size() const;
};

class DimensionIterator {
public:
	int index;
	Dimension* dim;
	/* return true if dimension index started again */
	bool next(Parameters* pars = NULL);
	DimensionIterator();
	double canPeek(int offset) const;
	double peek(int offset) const;
};

class SpaceIterator;

class ParamSpace {

public:
	std::vector<Dimension*> dimensions;

	SpaceIterator* current;

	size_t size() const;

	void add(Dimension* dim);

	size_t dim2ix(const SpaceIterator& si) const;

	void ix2dim(int ix, SpaceIterator& si) const;

	void iterate(Parameters& pars,std::function<void(const SpaceIterator&)> body);
	void iterate(std::function<void(const SpaceIterator&)> body);

	size_t currentInnerDimBase() const;

	ParamSpace() :current(NULL){};
};

class ParamSpaceValues {
public:
	ParamSpace* ps;

	Vector values;

	ParamSpaceValues(ParamSpace& ps);

	void initialize();

	double set(const SpaceIterator& si, double v);
	double get(const SpaceIterator& si) const;

	std::function<double(double)> dimInterpolator(int dimIx);

	void debug_dump();
};

class SpaceIterator {
public:
	bool next(Parameters* pars = NULL);

	std::vector<DimensionIterator> its;

	ParamSpace* ps;
	Parameters& par;

	SpaceIterator(ParamSpace* ps, Parameters& par);

};

//class DimInterpolate : public Fun {
//public:
//
//	const Dimension& dim;
//	const ParamSpaceValues& values;
//
//	DimInterpolate(const Dimension& dim, const ParamSpaceValues& values);
//	virtual double call(double v);
//
//};