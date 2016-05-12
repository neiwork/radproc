#include "SpaceIterator.h"

#include "ParamSpace.h"
#include "DimensionIterator.h"
#include "Dimension.h"

#include <algorithm>


SpaceIterator::SpaceIterator(const ParamSpace& ps, Parameters& par, std::initializer_list<int> fixedDimensions) :ps(ps), par(par), its(), coord(ps)
{
	std::vector<int> fixed = fixedDimensions;
	for (size_t i = 0; i < ps.dimensions.size(); ++i){
		int fix = i < fixed.size() ? fixed[i] : -1;
		its.push_back(DimensionIterator(ps.dimensions[i], coord.dims[i], fix));
		ps.dimensions[i]->update(par, std::max(0, fix));
	}
	its.front().reset();
}

bool SpaceIterator::next(Parameters* pars /*= NULL*/)
{
	size_t ii = 0;
	while (ii < its.size() && !its[ii].next(pars)) {
		++ii;
	}
	return ii < its.size();
}

SpaceCoord SpaceIterator::moved(std::initializer_list<int> deltas) const
{
	SpaceCoord newCoord = coord;
	size_t coord = 0;
	for (int delta : deltas) {
		newCoord[coord] += delta;
		coord++;
	}
	return newCoord;
}

bool SpaceIterator::canPeek(std::initializer_list<int> deltas) const
{
	size_t dim = 0;
	for (int delta : deltas) {
		if (!its[dim].canPeek(delta)) {
			return false;
		}
		dim++;
	}
	return true;
}


bool SpaceIterator::canPeek(const SpaceCoord& c) const
{
	for (size_t i = 0; i < its.size(); i++) {
		if (!its[i].canPeekAbs(c.dims[i])) {
			return false;
		}
	}
	return true;
}


