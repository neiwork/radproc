#include "SpaceCoord.h"

#include "ParamSpace.h"
#include "DimensionCoord.h"

SpaceCoord::SpaceCoord(const ParamSpace& ps)
{
	for (auto d:ps.dimensions){
		dims.push_back(0);
	}
}

SpaceCoord::SpaceCoord(std::initializer_list<DimensionCoord> dims):dims(dims)
{
}

SpaceCoord::SpaceCoord(const size_t dimensions)
{
	for (size_t i = 0; i < dimensions; ++i){
		dims.push_back(0);
	}
}

const DimensionCoord& SpaceCoord::operator[](const size_t& index) const
{
	return dims[index];
}

DimensionCoord& SpaceCoord::operator[](const size_t& index)
{
	return dims[index];
}


