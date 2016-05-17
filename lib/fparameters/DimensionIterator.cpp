#include "DimensionIterator.h"

#include "DimensionCoord.h"
#include "Dimension.h"

#include <algorithm>

DimensionIterator::DimensionIterator(Dimension* dim, DimensionCoord& coord, int fixed) :dim(dim), index(coord), fixed(fixed)
{
	if (fixed >= 0) {
		coord = fixed;
	}
}

bool DimensionIterator::next()
{
	index++;
	bool finished = fixed>=0?(index>fixed):(index == dim->size());
	if (finished) {
		// start over
		index = std::max(0,fixed);
	}
	return !finished;
}

double DimensionIterator::peek(int offset) const
{
	return dim->values[index + offset];
}

void DimensionIterator::reset()
{
	index = std::max(fixed,0)-1;
}

double DimensionIterator::canPeekAbs(int index) const
{
	return (0 <= index ) && (index  < (int)dim->size());
}

double DimensionIterator::canPeek(int offset) const
{
	return canPeekAbs(index + offset);
}


