#pragma once

#include "SpaceCoord.h"
#include "DimensionIterator.h"

class SpaceIterator {
public:
	bool next();

	SpaceCoord coord;

	std::vector<DimensionIterator> its;

	const ParamSpace& ps;

	SpaceIterator(const ParamSpace& ps, std::initializer_list<int> fixedDimensions = {});

	SpaceCoord moved(std::initializer_list<int> offsets) const;

	bool canPeek(std::initializer_list<int> offsets) const;
	bool canPeek(const SpaceCoord&) const;

	operator const SpaceCoord&() const {
		return coord;
	}

	double val(DimensionCoord dim) const;
};


