#include "vectorOperations.h"


#include <vector>
#include <algorithm>

double getMinValueVector(const Vector& v)  
{

    double minValue = *min_element (v.begin(), v.end());;

	return minValue;
}


double getMaxValueVector(const Vector& v)  
{

    double maxValue = *max_element (v.begin(), v.end());

	return maxValue;
}