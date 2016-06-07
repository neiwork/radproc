#pragma once


#include "DimensionCoord.h"
#include <fmath/interpolation.h>

#include <functional>

class SpaceCoord;
class ParamSpace;
class SpaceIterator;

class ParamSpaceValues {
public:
	const ParamSpace& ps;

	typedef std::function<void(std::ostream&, double)> DumpFun;

	Vector values;

	ParamSpaceValues(const ParamSpace& ps, bool initialize = true);
	ParamSpaceValues(const ParamSpace& ps, double fillValue );
	ParamSpaceValues(const ParamSpace& ps, std::function<double(const SpaceIterator& i)> initializer);
	ParamSpaceValues(const ParamSpaceValues& ps);

	ParamSpaceValues& operator =(const ParamSpaceValues& ps);

	void initialize();
	void initialize( double fillValue );

	void fill(const std::function<double(const SpaceIterator& i)>& f, std::initializer_list<int> fixed = {});

	double set(const SpaceCoord& si, double v);
	double get(const SpaceCoord& si) const;

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


