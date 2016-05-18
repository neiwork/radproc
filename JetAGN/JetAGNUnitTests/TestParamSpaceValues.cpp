#include "stdafx.h"
#include "CppUnitTest.h"

#include <fparameters/ParamSpaceValues.h>
#include <fparameters/ParamSpace.h>
#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>

#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace JetAGNUnitTests
{		

	void zeroToN(Vector& v) {
		for (int i = 0; i < (int)v.size(); ++i) {
			v[i] = i;
		}
	}

	TEST_CLASS(TestParamSpaceValues)
	{
	public:
		
		TEST_METHOD(MultiDimensionalInterpolate)
		{
			ParamSpace ps;
			ps.add(new Dimension(2, zeroToN));
			ps.add(new Dimension(3, zeroToN));
			ps.add(new Dimension(2, zeroToN));

			// create the PSVs
			ParamSpaceValues m(ps); // some PSV

			// iterate the dimensions, and for each combination of values; 
			// compute a value and store it appropriate location in PSV
			int c = 0;
			m.fill([&c](const SpaceIterator& i){
				return c++;
			});

			SpaceCoord p({ 0, 0, 0 });

			const double resultWithoutFallback = m.interpolate({ { 0, 0.1 }, { 1, 0.5 }, { 2, 0 } });
			const double resultWithFallback = m.interpolate({ { 0, 0.1 }, { 1, 0.5 } }, &p);
			const double correctResult{ 1.0999999999999999 };
			Assert::AreEqual(correctResult, resultWithoutFallback, L"interpolate did not work correctly without fallback coordinates.", LINE_INFO());
			Assert::AreEqual(correctResult, resultWithFallback, L"interpolate did not work correctly using fallback coordinates.", LINE_INFO());
		}

	};
}