#include "ParamSpace.h"

#include "SpaceIterator.h"
#include "Dimension.h"
#include "SpaceCoord.h"

#include <omp.h>

void ParamSpace::iterate(std::function<void(const SpaceIterator&)> body, std::initializer_list<int> fixedDimensions) const
{
	SpaceIterator it(*this,fixedDimensions);
	while (it.next()) {
		body(it);
	}
}

void ParamSpace::parallelize(std::function<void(const SpaceIterator&)> body) const
{
	const size_t sz{ size() };
	#pragma omp parallel
	for (size_t i = 0; i < sz; ++i) {
		SpaceIterator it(*this);
		ix2dim(i, it.coord);
		body(it);
	}
}

void ParamSpace::add(Dimension* dim)
{
	dimensions.push_back(dim);
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
		ix /= dimensions[i]->size();
	}
}

ParamSpace::ParamSpace()
{

}

const Dimension& ParamSpace::operator[](const size_t& index) const
{
	return *dimensions[index];
}
