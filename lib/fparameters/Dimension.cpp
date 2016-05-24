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

const double& Dimension::operator[](const int& index) const
{
	return values[index];
}


