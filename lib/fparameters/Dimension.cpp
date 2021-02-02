#include "Dimension.h"

#include "ParamSpaceValues.h"
#include "ParamSpace.h"

#include <iostream>

int Dimension::size() const
{
	return values.size();
}

//void Dimension::update(Parameters& pars, size_t index)
//{
//	// update the parameters object using the pointer-to-member
//	//if (fieldRef != nullptr) {
//	//	pars.*fieldRef = values[index];
//	//}
//}

Dimension::Dimension(int size, std::function<void(Vector&)> initializer)
{
	values.resize(size, 0.0);
	initializer(values);
}

Dimension::Dimension(int size, std::function<void(Vector&)> initializer, std::function<double(double)> to_linear, std::function<double(double)> from_linear) 
	:to_linear(to_linear), from_linear(from_linear)
{
	values.resize(size, 0.0);
	initializer(values);
}

const double& Dimension::operator[](const int& index) const
{
	return values[index];
}

double Dimension::interpolateIndex(double v, int l, int u) const {
	if(u == l) {
		return l;
	}
	double value = to_linear(v);
	double lower = to_linear(values[l]);
	double upper = to_linear(values[u]);
	return (value - lower) / (upper - lower) + l;
}
