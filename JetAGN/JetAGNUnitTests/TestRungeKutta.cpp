#include "CppUnitTest.h"

#include <fmath\RungeKutta.h>



using namespace Microsoft::VisualStudio::CppUnitTestFramework;
//using boost::property_tree::ptree;

namespace RungeKuttaTests
{

	double min(double x)
	{
		return x;
	}
	double max(double x)
	{
		return 2.0*x;
	}
	double f(double x, double y)
	{
		return 3 * x;
	}

	double tol = 5.0e-1;

	TEST_CLASS(Integrators)
	{
	public:

		TEST_METHOD(TestDoubleIntegral)
		{
			double integral = RungeKutta(1.0, 2.0, min, max, f);
			double expected = 7.0;
			

			//Assert.That(expected, Is.EqualTo(actual).Within(tolerance);
			Assert::AreEqual(expected, integral, tol, L"Double integration is not working", LINE_INFO());
		}

		TEST_METHOD(TestSingleIntegral)
		{
			double integral = RungeKuttaSimple(2.0, 4.0, max);
			double expected = 12.0;

			Assert::AreEqual(expected, integral, tol, L"Simple integration is not working", LINE_INFO());
			
	
		}

	};
}