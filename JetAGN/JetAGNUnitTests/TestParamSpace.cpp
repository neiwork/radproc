#include "CppUnitTest.h"

#include <fparameters/ParamSpaceValues.h>
#include <fparameters/ParamSpace.h>
#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>
#include <JetAGN/checks.h>

#include <numeric>
#include <set>
#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace std;
namespace JetAGNUnitTests
{

	TEST_CLASS(TestParamSpace)
	{
	public:

		TEST_METHOD(TestIxConversion)
		{
			ParamSpace ps;
			ps.add(new Dimension(2, zeroToN));
			ps.add(new Dimension(3, zeroToN));
			ps.add(new Dimension(4, zeroToN));

			const size_t sz{ ps.size() };
			typedef tuple<int, int, int> T3;
			set<T3> visited;
			SpaceCoord c(ps);
			for (size_t i = 0; i < sz; ++i) {
				ps.ix2dim(i, c);
				Assert::AreEqual(i, ps.dim2ix(c), L"Inverse should work", LINE_INFO());
				T3 t{ c[0], c[1], c[2] };
				Assert::IsFalse(visited.count(t)>0,L"Coordinate already mapped",LINE_INFO());
				visited.insert(t);
			}
			Assert::AreEqual(ps.size(),visited.size(), L"Should have all coordinates", LINE_INFO());
		}
	};
}