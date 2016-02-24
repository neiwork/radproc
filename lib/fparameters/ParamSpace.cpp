#include "ParamSpace.h"
#include <fmath/interpolation.h>
#include <iterator>
#include <iostream>
#include <iomanip>

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
	if (par != &Parameters::E || psv.ps.dimensions[0] != this) {
		std::cout << "Por ahora solo se puede interpolar la dimension 0, que debe ser la energia." << std::endl;
		throw;
	}

	const Vector& key = values;
	const Vector& val = psv.values;

	size_t base = psv.ps.currentInnerDimBase();
	return interpol(x, values, psv.values, values.size(), base);
}


DimensionIterator::DimensionIterator(Dimension* dim, DimensionCoord& coord, int fixed) :dim(dim), index(coord), fixed(fixed)
{
	if (fixed >= 0) {
		coord = fixed;
	}
}

bool DimensionIterator::next(Parameters* pars /*= NULL*/)
{
	index++;
	bool finished = fixed>=0?(index>fixed):(index == dim->size());
	if (finished) {
		// start over
		index = std::max(0,fixed);
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

SpaceCoord::SpaceCoord(const ParamSpace& ps)
{
	for (size_t i = 0; i < ps.dimensions.size(); ++i){
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

void ParamSpace::updateDerivations(const SpaceIterator& i) const {
	for (auto d : derivations) {
		d(i);
	}
}

void ParamSpace::iterate(Parameters& p, std::function<void(const SpaceIterator&)> body, std::initializer_list<int> fixedDimensions) const
{
	SpaceIterator it(*this,p,fixedDimensions);
	current = &it;
	while (it.next(&p)) {
		updateDerivations(it);
		body(it);
	}
	current = NULL;
}

void ParamSpace::iterate(std::function<void(const SpaceIterator&)> body, std::initializer_list<int> fixedDimensions) const
{
	iterate(Parameters(parameters), body, fixedDimensions);
}

void ParamSpace::add(Dimension* dim)
{
	dimensions.push_back(dim);
}

void ParamSpace::addDerivation( const DerivationFun& f )
{
	derivations.push_back(f);
}

size_t ParamSpace::size() const
{
	size_t N = 1;
	for (auto v : dimensions) { N *= v->size(); }
	return N;
}

size_t ParamSpace::dim2ix(const SpaceCoord& coord) const
{
	size_t ix = coord[0];
	size_t o = dimensions[0]->size();
	for (size_t i = 1; i < coord.dims.size(); ++i) {
		ix += coord[i] * o;
		o *= dimensions[i]->size();
	}
	return ix;
}

void ParamSpace::ix2dim(int ix, SpaceCoord& si) const
{
	for (size_t i = 0; i < dimensions.size(); ++i) {
		si[i] = ix % dimensions[i]->size();
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
	values.resize(ps.size());
}

ParamSpaceValues::ParamSpaceValues(const ParamSpace& ps) :ps(ps)
{
}


ParamSpaceValues::ParamSpaceValues(const ParamSpace& ps, std::function<double(const SpaceIterator& i)> initializer) : ps(ps)
{
	initialize();
	fill(initializer);
}

double ParamSpaceValues::interpolate(std::initializer_list<double> dimValues) const {
	std::vector<double> values(dimValues), values_idx;
	const size_t D = ps.dimensions.size();
	if (dimValues.size() != D) {
		std::cout << "ERROR: Invalid number of coordinate for this ParamSpaceValues." << std::endl;
		throw;
	}
	
	// Calculo upper y lower que son las esquinas del cubo alrededor del punto a interpolar
	std::vector<size_t> lower, upper;
	for (size_t i = 0; i < D; ++i) {
		auto& dimv = ps[i].values;
		{	
			auto v = values[i];

			auto lb = std::lower_bound(dimv.begin(), dimv.end(), v);
			// {  { *(ub-1) < v && v <= *ub }

			size_t u = lb - dimv.begin(); // truco de C++ convierte de std::iterator a indice
			if (lb != dimv.end()) {
				int l;
				if (*lb != v) {
					l = u - 1;
				} else {
					l = u;
				}
				lower.push_back(l);
				upper.push_back(u);
				
				if (u != l) {
					values_idx.push_back((v - dimv[l]) / (dimv[u] - dimv[l]) + l); // values in 'index space'
				} else {
					values_idx.push_back(l); // values in 'index space'
				}
			}
			else {
				std::cout << " out of range" << std::endl;
				throw;
			}
		}
	}
	
	// DEBUG: print upper and lower bounds (indexes)
/*	for (size_t i = 0; i < D; ++i) {
		std::cout << lower[i] << " < ";
		std::cout << values_idx[i] << " < ";
		std::cout << upper[i] << std::endl;
	}*/
	double value = 0.0;
	
	// iterate all possible hyperquadrants defined by values, compute area(weight) and add to value
	double dbg_total_area = 0;
	for (size_t i = 0; i < std::pow(2, D); ++i) {
		SpaceCoord corner(ps), opposite(ps);
//		std::cout << "quadrant ";
		size_t quads = i;
		for (size_t j = 0; j < D; ++j, quads /= 2) {
//			std::cout << quads % 2;
			corner[j] = quads % 2 ? lower[j] : upper[j];
			opposite[j] = quads % 2 ? upper[j] : lower[j];
		}
		// {corner now is the space coordinate corresponding to the current corner of the hypercube}

		// now, compute the area of the hyperquadrant opposite to this corner
		double area = 1;
		for (size_t j = 0; j < D; j++) {
			if (opposite[j] != corner[j]) {
				area *= std::abs(values_idx[j] - opposite[j]);
			} else {
				area /= 2;
			}
		}
		dbg_total_area += area;
		value += area * get(corner);
//		std::cout << " w:" << area << ", v:" << get(corner) << std::endl;
	}

//	std::cout << "dbgarea:" << dbg_total_area << std::endl;
//	std::cout << "value:" << value << std::endl;

	return value;
}


double ParamSpaceValues::get(const SpaceCoord& si) const
{
		return values[ps.dim2ix(si)];
}

double ParamSpaceValues::set(const SpaceCoord& si, double v)
{
	return values[ps.dim2ix(si)] = v;
}

void ParamSpaceValues::fill(const std::function<double(const SpaceIterator& i)>& f)
{
	ps.iterate([&f, this](const SpaceIterator& i){
		this->set(i, f(i));
	});
}

void ParamSpaceValues::dump(std::ostream& o, DumpFun dumper)
{
	for (size_t i = 0; i < values.size(); ++i) {
		dumper(o,values[i]);
		if ((i + 1) % (ps.dimensions.front()->size()) == 0)
			std::cout << std::endl;
		else
			std::cout << " ";
	}
}

void ParamSpaceValues::std_dump(std::ostream& o, double v){
	o << std::setw(2) << v;
};

void ParamSpaceValues::debug_dump()
{
	dump(std::cout, std_dump);
}
std::function<double(double)> ParamSpaceValues::dimInterpolator(int dimIx)
{
	Dimension* d = ps.dimensions[dimIx];
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

SpaceIterator* ParamSpace::current = NULL;
