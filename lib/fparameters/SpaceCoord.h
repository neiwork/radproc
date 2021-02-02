#pragma once

#include "DimensionCoord.h"
#include <vector>

#include <cstddef>

class ParamSpace;

class SpaceCoord {
public:
	std::vector<DimensionCoord> dims;        // dims is an int vector

    // Constructors
	SpaceCoord(const ParamSpace& ps);    // ps store the address of a ParamSpace type
	explicit SpaceCoord(const size_t dimensions);     // size_t is an integer type used for indexes
	SpaceCoord(std::initializer_list<DimensionCoord> dims);

	const DimensionCoord& operator[](const size_t& index) const;  // operator[]() store the address of an int

	DimensionCoord& operator[](const size_t& index);        // operator[]() store the address of an int
};

