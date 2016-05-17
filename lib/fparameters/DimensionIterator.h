#pragma once

#include "DimensionCoord.h"

class Dimension;

class DimensionIterator {
public:
	DimensionCoord &index;
	Dimension* dim;
	int fixed;
	/* return true if dimension index started again */
	bool next();
	DimensionIterator(Dimension*,DimensionCoord&,int fixed=-1);
	double canPeek(int offset) const;
	double canPeekAbs(int index) const;
	double peek(int offset) const;
	void reset();
	inline double val() const {
		return peek(0);
	}
};


