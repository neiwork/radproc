#include "ParamSpace.h"
#include <fmath/interpolation.h>
#include <iterator>
#include <iostream>

size_t Dimension::size() const
{
	return values.size();
}

void Dimension::update(Parameters& pars, size_t index)
{
	// update the parameters object using the pointer-to-member
	pars.*par = values[index];
}

Dimension::Dimension(size_t size, double Parameters::*par, std::function<void(Vector&)> initializer) :par(par)
{
	values.resize(size, 0.0);
	initializer(values);
}

double Dimension::interpolate(double x, const ParamSpaceValues& psv) const
{
	if (par != &Parameters::E || psv.ps->dimensions[0] != this) {
		std::cout << "Por ahora solo se puede interpolar la dimension 0, que debe ser la energia." << std::endl;
		throw;
	}

	const Vector& key = values;
	const Vector& val = psv.values;

	size_t base = psv.ps->currentInnerDimBase();
	return interpol(x, values, psv.values, values.size(), base);
}


DimensionIterator::DimensionIterator() :dim(NULL), index(-1)
{

}

bool DimensionIterator::next(Parameters* pars /*= NULL*/)
{
	index++;
	bool finished = (index == dim->size());
	if (finished) {
		// start over
		index = 0;
	}
	if (pars) {
		dim->update(*pars, index);
	}
	return !finished;
}

double DimensionIterator::peek(int offset) const
{
	return dim->values[index + offset];
}

double DimensionIterator::canPeek(int offset) const
{
	return (0 <= index + offset) && (index + offset < (int)dim->size());
}

SpaceIterator::SpaceIterator(ParamSpace* ps, Parameters& par) :ps(ps),par(par),its(ps->dimensions.size())
{
	for (size_t i = 0; i < ps->dimensions.size(); ++i){
		its[i].dim = ps->dimensions[i];
		its[i].index = 0;
		ps->dimensions[i]->update(par, 0);
	}
	its.front().index = -1;
}

bool SpaceIterator::next(Parameters* pars /*= NULL*/)
{
	size_t ii = 0;
	while (ii < its.size() && !its[ii].next(pars)) {
		++ii;
	}
	return ii < its.size();
}

void ParamSpace::iterate(Parameters& p, std::function<void(const SpaceIterator&)> body)
{
	SpaceIterator it(this,p);
	current = &it;
	while (it.next(&p)) {
		body(it);
	}
	current = NULL;
}

void ParamSpace::iterate(std::function<void(const SpaceIterator&)> body)
{
	iterate(Parameters(parameters), body);
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

size_t ParamSpace::dim2ix(const SpaceIterator& si) const
{
	size_t ix = si.its.front().index;
	for (size_t i = 1; i < si.its.size(); ++i) {
		ix += si.its[i].index * si.its[i-1].dim->size();
	}
	return ix;
}

void ParamSpace::ix2dim(int ix, SpaceIterator& si) const
{
	for (auto i = si.its.begin(); i != si.its.end(); ++i) {
		i->index = ix % i->dim->size();
	}
}

size_t ParamSpace::currentInnerDimBase() const
{
	// find the current location; and then substract the innermost dimension index.
	// this is the starting point for energy-based interpolation of PSVs
	if (current && dimensions.size() > 1) {
		return dim2ix(*current) - current->its[0].index;
	} else {
		return 0;
	}
}

void ParamSpaceValues::initialize()
{
	values.resize(ps->size());
}

ParamSpaceValues::ParamSpaceValues(ParamSpace& ps) :ps(&ps)
{
}

double ParamSpaceValues::get(const SpaceIterator& si) const
{
	return values[ps->dim2ix(si)];
}

double ParamSpaceValues::set(const SpaceIterator& si, double v)
{
	return values[ps->dim2ix(si)] = v;
}

void ParamSpaceValues::debug_dump()
{
	for (size_t i = 0; i < values.size(); ++i) {
		std::cout << values[i];
		if ((i + 1) % (ps->dimensions.front()->size()) == 0)
			std::cout << std::endl;
		else
			std::cout << " ";
	}
}

std::function<double(double)> ParamSpaceValues::dimInterpolator(int dimIx)
{
	Dimension* d = ps->dimensions[dimIx];
	ParamSpaceValues* psv = this;
	return [d, psv](double x) {
		return d->interpolate(x, *psv);
	};
}

//double DimInterpolate::call(double v) {
//	return dim.interpolate(v, values);
//}
//
//DimInterpolate::DimInterpolate(const Dimension& dim, const ParamSpaceValues& values) :dim(dim), values(values)
//{
//
//}
