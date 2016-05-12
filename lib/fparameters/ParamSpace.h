#pragma once
#include "parameters.h"

class SpaceCoord;
class Dimension;
class SpaceIterator;

class ParamSpace {

	void updateDerivations(const SpaceIterator& i) const;

public:

	typedef std::function<void(const SpaceIterator&)> DerivationFun;

	std::vector<Dimension*> dimensions;
	std::vector<const DerivationFun> derivations;

	size_t size() const;

	void add(Dimension* dim);

	void addDerivation( const DerivationFun& f );

	size_t dim2ix(const SpaceCoord& si) const;

	void ix2dim(int ix, SpaceCoord& si) const;

	void iterate(Parameters& pars, std::function<void(const SpaceIterator&)> body, std::initializer_list<int> fixedDimensions = {}) const;
	void iterate(std::function<void(const SpaceIterator&)> body, std::initializer_list<int> fixedDimensions = {}) const;

	ParamSpace(){};

	const Dimension& operator[](const size_t& index) const {
		return *dimensions[index];
	}
};

