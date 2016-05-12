#pragma once

#include "SpaceCoord.h"
#include "DimensionIterator.h"

class SpaceIterator {
public:
	bool next(Parameters* pars = nullptr);

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


