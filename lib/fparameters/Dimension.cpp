#include "Dimension.h"

#include "ParamSpaceValues.h"
#include "ParamSpace.h"

#include <iostream>

size_t Dimension::size() const
{
	return values.size();
}

void Dimension::update(Parameters& pars, size_t index)
{
	// update the parameters object using the pointer-to-member
	if (fieldRef != nullptr) {
		pars.*fieldRef = values[index];
	}
}

Dimension::Dimension(size_t size, std::function<void(Vector&)> initializer, double Parameters::*fieldRef) :fieldRef(fieldRef)
{
	values.resize(size, 0.0);
	initializer(values);
}

double Dimension::interpolate(double x, const ParamSpaceValues& psv) const
{
	if (fieldRef != &Parameters::E || psv.ps.dimensions[0] != this) {
		std::cout << "Por ahora solo se puede interpolar la dimension 0, que debe ser la energia." << std::endl;
		throw;
	}

	const Vector& key = values;
	const Vector& val = psv.values;

	size_t base = psv.ps.currentInnerDimBase();
	return interpol(x, values, psv.values, values.size(), base);
}

const double& Dimension::operator[](const size_t& index) const
{
	return values[index];
}


