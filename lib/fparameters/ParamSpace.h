#pragma once
#include "parameters.h"

#include <fmath/mathematics.h>

#include <vector>
#include <functional>
#include <numeric>
#include <algorithm>
#include <ostream>

typedef int DimensionCoord;

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

class DimensionIterator {
public:
	DimensionCoord &index;
	Dimension* dim;
	int fixed;
	/* return true if dimension index started again */
	bool next(Parameters* pars = NULL);
	DimensionIterator(Dimension*,DimensionCoord&,int fixed=-1);
	double canPeek(int offset) const;
	double canPeekAbs(int index) const;
	double peek(int offset) const;
	void reset();
};

class SpaceIterator;
class ParamSpace;

class SpaceCoord {
public:
	std::vector<DimensionCoord> dims;

	SpaceCoord(const ParamSpace& ps);
	SpaceCoord(const size_t dimensions);
	SpaceCoord(std::initializer_list<DimensionCoord> dims);

	const DimensionCoord& operator[](const size_t& index) const;

	DimensionCoord& operator[](const size_t& index);
};

class ParamSpace {

	void updateDerivations(const SpaceIterator& i) const;

public:

	typedef std::function<void(const SpaceIterator&)> DerivationFun;

	std::vector<Dimension*> dimensions;
	std::vector<const DerivationFun> derivations;

	static SpaceIterator* current;

	size_t size() const;

	void add(Dimension* dim);

	void addDerivation( const DerivationFun& f );


	size_t dim2ix(const SpaceCoord& si) const;

	void ix2dim(int ix, SpaceCoord& si) const;

	void iterate(Parameters& pars, std::function<void(const SpaceIterator&)> body, std::initializer_list<int> fixedDimensions = {}) const;
	void iterate(std::function<void(const SpaceIterator&)> body, std::initializer_list<int> fixedDimensions = {}) const;

	size_t currentInnerDimBase() const;

	ParamSpace(){};

	const Dimension& operator[](const size_t& index) const {
		return *dimensions[index];
	}
};

class ParamSpaceValues {
public:
	const ParamSpace& ps;

	typedef std::function<void(std::ostream&, double)> DumpFun;

	Vector values;

	ParamSpaceValues(const ParamSpace& ps, bool initialize = true);
	ParamSpaceValues(const ParamSpace& ps, std::function<double(const SpaceIterator& i)> initializer);
	ParamSpaceValues(const ParamSpaceValues& ps);

	ParamSpaceValues& operator =(const ParamSpaceValues& ps);

	void initialize();

	void fill(const std::function<double(const SpaceIterator& i)>& f, std::initializer_list<int> fixed = {});

	double set(const SpaceCoord& si, double v);
	double get(const SpaceCoord& si) const;

	std::function<double(double)> dimInterpolator(int dimIx);

	class InterpolateDim {
	public:
		DimensionCoord dim;
		double value;
		InterpolateDim(DimensionCoord d, double v) :dim(d), value(v) {}
	};

	double interpolate(std::initializer_list<InterpolateDim> dimValues, const SpaceCoord* fallback=0) const;

	// debug:
		void debug_dump();
		void dump(std::ostream& o, DumpFun dumper);
		static void std_dump(std::ostream& o, double v);
};


class SpaceIterator {
public:
	bool next(Parameters* pars = NULL);

	SpaceCoord coord;

	std::vector<DimensionIterator> its;

	const ParamSpace& ps;
	Parameters& par;

	SpaceIterator(const ParamSpace& ps, Parameters& par, std::initializer_list<int> fixedDimensions = {});

	SpaceCoord moved(std::initializer_list<int> offsets) const;

	bool canPeek(std::initializer_list<int> offsets) const;
	bool canPeek(const SpaceCoord&) const;

	operator const SpaceCoord&() const {
		return coord;
	}

};
