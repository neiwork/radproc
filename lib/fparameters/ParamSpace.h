#pragma once

#include <vector>
#include <functional>

class SpaceCoord;
class Dimension;
class SpaceIterator;

class ParamSpace {
public:

	std::vector<Dimension*> dimensions;

	size_t size() const;

	void add(Dimension* dim);

	size_t dim2ix(const SpaceCoord& si) const;

	void ix2dim(int ix, SpaceCoord& si) const;

	void iterate(std::function<void(const SpaceIterator&)> body, std::initializer_list<int> fixedDimensions = {}) const;
	void parallelize(std::function<void(const SpaceIterator&)> body) const;

	ParamSpace();

	const Dimension& operator[](const size_t& index) const;
};

