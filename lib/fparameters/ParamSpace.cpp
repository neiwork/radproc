#include "ParamSpace.h"

#include "SpaceIterator.h"
#include "Dimension.h"
#include "SpaceCoord.h"

//#include <fmath/interpolation.h>
//#include <iterator>
//#include <iostream>
//#include <iomanip>
//#include <stdexcept>
//#include <sstream>


void ParamSpace::updateDerivations(const SpaceIterator& i) const {
	for (auto d : derivations) {
		d(i);
	}
}

void ParamSpace::iterate(Parameters& p, std::function<void(const SpaceIterator&)> body, std::initializer_list<int> fixedDimensions) const
{
	SpaceIterator it(*this,p,fixedDimensions);
	current = &it;
	while (it.next(&p)) {
		updateDerivations(it);
		body(it);
	}
	current = NULL;
}

void ParamSpace::iterate(std::function<void(const SpaceIterator&)> body, std::initializer_list<int> fixedDimensions) const
{
	iterate(Parameters(parameters), body, fixedDimensions);
}

void ParamSpace::add(Dimension* dim)
{
	dimensions.push_back(dim);
}

void ParamSpace::addDerivation( const DerivationFun& f )
{
	derivations.push_back(f);
}

size_t ParamSpace::size() const
{
	size_t N = 1;
	for (auto v : dimensions) { N *= v->size(); }
	return N;
}

size_t ParamSpace::dim2ix(const SpaceCoord& coord) const
{
	size_t ix = coord[0];
	size_t o = dimensions[0]->size();
	for (size_t i = 1; i < coord.dims.size(); ++i) {
		ix += coord[i] * o;
		o *= dimensions[i]->size();
	}
	return ix;
}

void ParamSpace::ix2dim(int ix, SpaceCoord& si) const
{
	for (size_t i = 0; i < dimensions.size(); ++i) {
		si[i] = ix % dimensions[i]->size();
	}
}

size_t ParamSpace::currentInnerDimBase() const
{
	// find the current location; and then substract the innermost dimension index.
	// this is the starting point for energy-based interpolation of PSVs
	if (current && dimensions.size() > 1) {
		return dim2ix(*current) - current->its[0].index;
	} else {
		return 0;
	}
}


//double DimInterpolate::call(double v) {
//	return dim.interpolate(v, values);
//}
//
//DimInterpolate::DimInterpolate(const Dimension& dim, const ParamSpaceValues& values) :dim(dim), values(values)
//{
//
//}

SpaceIterator* ParamSpace::current = NULL;
