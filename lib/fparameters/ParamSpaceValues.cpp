#include "ParamSpaceValues.h"

#include "ParamSpace.h"
#include "DimensionCoord.h"
#include "SpaceCoord.h"
#include "Dimension.h"
#include "SpaceIterator.h"

#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>

void ParamSpaceValues::initialize()
{
	values.resize(ps.size());
}

void ParamSpaceValues::initialize(double fillValue)
{
	values.resize(ps.size(),fillValue);
}

ParamSpaceValues::ParamSpaceValues(const ParamSpace& ps, bool initialize):ps(ps)
{
	if (initialize) {
		this->initialize();
	}
}

ParamSpaceValues::ParamSpaceValues(const ParamSpace& ps, double fillValue) :ps(ps)
{
	initialize(fillValue);
}


ParamSpaceValues::ParamSpaceValues(const ParamSpace& ps, std::function<double(const SpaceIterator& i)> initializer) : ps(ps)
{
	initialize();
	fill(initializer);
}

double ParamSpaceValues::interpolate(std::initializer_list<InterpolateDim> dimValues, const SpaceCoord* fallback) const {
	const DimensionCoord D = ps.dimensions.size();

	if (!fallback && dimValues.size() != D) {
		std::ostringstream errormsg;
		errormsg << "Invalid number of coordinates " << dimValues.size() << "  for interpolating in " << D << " dimensions." << std::endl;
		throw std::runtime_error(errormsg.str());
	}

	std::vector<DimensionCoord> lower(D), upper(D);
	std::vector<double> values(D); // values
	std::vector<double> values_idx(D); // values in 'index space'

	if (fallback) { // if we have a fallback coord, fill in fixed && equal (no interpolation) coords.
		lower = fallback->dims;
		upper = fallback->dims;
		for (DimensionCoord d = 0; d < D; ++d) {
			values_idx[d] = lower[d];
			values[d] = ps[d][lower[d]];
		}
	}
	
	// Calculo upper y lower que son las esquinas del cubo alrededor del punto a interpolar
	for (InterpolateDim id:dimValues) {
		auto d = id.dim;
		auto& dimv = ps[d].values; // get the array of values defined for this particular Dimension
		{	
			auto v = values[d] = id.value;

			auto lb = std::lower_bound(dimv.begin(), dimv.end(), v);
			// {  { *(ub-1) < v && v <= *ub }

			size_t u = lb - dimv.begin(); // truco de C++ convierte de std::iterator a indice
			if (lb != dimv.end()) {
				int l;
				if (*lb != v) {
					l = u - 1;
				} else {
					l = u;
				}
				lower[d] = l;
				upper[d] = u;
			} else {
				std::ostringstream errormsg;
				errormsg << "Interpolation value " << id.value << " is not in range [ " << dimv.front() << " .. " << dimv.back() << " ] of dimension " << d << "." << std::endl;
 				throw std::runtime_error(errormsg.str());
			}
		}
	}

	// compute values for all dimensions in index space (possibly non-integers)

	for (DimensionCoord d = 0; d < D; ++d) {
		auto& dimv = ps[d].values; // get the array of values defined for this particular Dimension
		int l = lower[d];
		int u = upper[d];
		values_idx[d] = ps[d].interpolateIndex(values[d], l, u); // dim values in 'index space'
	}

	// DEBUG: print upper and lower bounds (indexes)
	//for (DimensionCoord i = 0; i < D; ++i) {
	//	std::cout << lower[i] << " < ";
	//	std::cout << values_idx[i] << " < ";
	//	std::cout << upper[i] << std::endl;
	//}
	double value = 0.0;
	
	const auto &to_linear = ps[(*(dimValues.begin())).dim].to_linear;
	const auto &from_linear = ps[(*(dimValues.begin())).dim].from_linear;
	
	// iterate all possible hyperquadrants defined by values, compute area(weight) and add to value
	double dbg_total_area = 0;
	for (unsigned int quadrant = 0; quadrant < std::pow(2, D); ++quadrant) {
		SpaceCoord corner(ps), opposite(ps);
		//std::cout << "quadrant ";
		// bits determine orientation in each dimension
		unsigned int orientation = quadrant; 
		for (DimensionCoord j = 0; j < D; ++j, orientation /= 2) {
			//std::cout << quads % 2;
			corner[j]   = orientation % 2 ? lower[j] : upper[j];
			opposite[j] = orientation % 2 ? upper[j] : lower[j];
		}
		// {corner now is the space coordinate corresponding to the current corner of the hypercube}

		// now, compute the area of the hyperquadrant opposite to this corner
		double area = 1;
		for (DimensionCoord j = 0; j < D; j++) {
			if (opposite[j] != corner[j]) {
				area *= std::abs(values_idx[j] - opposite[j]);
			} else {
				area /= 2;
			}
		}
		dbg_total_area += area;
		value += area * to_linear(get(corner));
		//std::cout << " w:" << area << ", v:" << get(corner) << std::endl;
	}

	//std::cout << "dbgarea:" << dbg_total_area << std::endl;
	//std::cout << "value:" << value << std::endl;

	return from_linear(value);
}

ParamSpaceValues::ParamSpaceValues(const ParamSpaceValues& psv):ps(psv.ps),values(psv.values) {
	
}

ParamSpaceValues& ParamSpaceValues::operator =(const ParamSpaceValues& psv) {
	// solo funciona bien si el PS actual es igual al PS del psv q se quiere asignar
	values = psv.values;
	return *this;
}


double ParamSpaceValues::get(const SpaceCoord& si) const
{
		return values[ps.dim2ix(si)];
}

double ParamSpaceValues::set(const SpaceCoord& si, double v)
{
	return values[ps.dim2ix(si)] = v;
}

void ParamSpaceValues::fill(const std::function<double(const SpaceIterator& i)>& f, std::initializer_list<int> fixed)
{
	ps.iterate([&f, this](const SpaceIterator& i){
		this->set(i, f(i));
	},fixed);
}

void ParamSpaceValues::dump(std::ostream& o, DumpFun dumper)
{
	for (size_t i = 0; i < values.size(); ++i) {
		dumper(o,values[i]);
		if ((i + 1) % (ps.dimensions.front()->size()) == 0)
			std::cout << std::endl;
		else
			std::cout << " ";
	}
}

void ParamSpaceValues::std_dump(std::ostream& o, double v){
	o << std::setw(2) << v;
};

void ParamSpaceValues::debug_dump()
{
	dump(std::cout, std_dump);
}