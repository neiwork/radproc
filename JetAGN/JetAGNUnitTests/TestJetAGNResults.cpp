#include "CppUnitTest.h"

#include <fparameters/Dimension.h>

#include <JetAGN/distribution.h>
#include <JetAGN/injection.h>
#include <JetAGN/State.h>
#include <JetAGN/modelParameters.h>
#include <JetAGN/checkPower.h>
#include <JetAGN/checks.h>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using boost::property_tree::ptree;

namespace JetAGNUnitTests
{
	/* Default parameter config for testing. */
	boost::property_tree::ptree config() {

		boost::property_tree::ptree cfg, dim;
		dim.add<int>("energy.samples", 4);
		dim.add<int>("radius.samples", 3);
		dim.add<int>("time.samples", 3);
		cfg.add_child("particle.default.dim", dim);

		cfg.add<double>("particle.electron.mass", electronMass);
		cfg.add<double>("particle.electron.dim.energy.min", 6.0);
		cfg.add<double>("particle.electron.dim.energy.max", 15.0);
		cfg.add<double>("particle.photon.dim.energy.min", -6.0);
		cfg.add<double>("particle.photon.dim.energy.max", 12.0);

		return std::move(cfg);
	}

	TEST_CLASS(JetAGNResults)
	{
	public:

		TEST_METHOD(TestStateConfiguration)
		{
			ptree cfg{ config() };
			setParameters(cfg);
			State model(cfg);
			
			Assert::AreEqual<size_t>(4, model.electron.ps[0].size(), L"energy dimension size", LINE_INFO());
			Assert::AreEqual<size_t>(3, model.electron.ps[1].size(), L"energy dimension size", LINE_INFO());
			Assert::AreEqual<size_t>(3, model.electron.ps[2].size(), L"energy dimension size", LINE_INFO());
		}

		TEST_METHOD(TestEnergyDimensionValues11)
		{
			ptree cfg{config()};
			cfg.put<int>("particle.default.dim.energy.samples", 11);

			setParameters(cfg);
			State model(cfg);
			
			// check construction of standard 11-value energy vector
			std::vector<double> expected{ { 1.6000000000000001e-006, 1.6000000000000003e-005, 0.00016000000000000007, 0.0016000000000000009, 0.016000000000000011, 0.16000000000000014, 1.6000000000000016, 16.000000000000018, 160.00000000000020, 1600.0000000000023, 16000.000000000025 } };

			Assert::IsTrue(check_vec(expected, model.electron.ps[0].values), L"wrong energy dimension vector", LINE_INFO());
		}


		TEST_METHOD(TestInjection)
		{
			ptree cfg{config()};
			setParameters(cfg);
			State model(cfg);

			injection(model.electron, model);

			double injectedPower = computeInjectedPower(model.electron.injection, 0);
			Assert::AreEqual(1.2898496718677250e+034, injectedPower, L"wrong injected power", LINE_INFO());

		}

		TEST_METHOD(TestDistribution)
		{
			// OVERRIDE DISCRETIZATION TO MAKE SURE TEST RUNS FAST
			// >>
			ptree cfg{config()};
			setParameters(cfg);
			State model(cfg);
			// <<

			injection(model.electron, model);
			distribution(model.electron, model);

			// HARDCODED EXPECTED RESULTS IN THE PSV (36 values)
			std::vector<double> expected{ { 2.7350568658788040e-007, 6.3622294948870793e-012, 2.4501820938773908e-023, 0.00000000000000000, 8.6490092262586863e-009, 2.0119136200545765e-013, 7.7481560988130572e-025, 0.00000000000000000, 2.7350568658788044e-010, 6.3622294948870794e-015, 2.4501820938773911e-026, 0.00000000000000000, 2.7351350281398751e-007, 6.3642955727248140e-012, 2.4683661547834183e-021, 0.00000000000000000, 8.6494266676343737e-009, 2.0135548889937072e-013, 2.2343163115867641e-021, 0.00000000000000000, 8.1161568079957346e-008, 2.5697611275525037e-012, 8.2407130680496509e-022, 0.00000000000000000, 2.7429311002179792e-007, 6.5703706212366638e-012, 2.4622465256150287e-019, 0.00000000000000000, 9.1649148042644592e-009, 4.9265672526681772e-013, 6.2257382988581396e-018, 0.00000000000000000, 4.4941289544580628e-010, 1.0999721507612061e-013, 2.5536564239863357e-017, 0.00000000000000000 } };

			Assert::IsTrue(check_vec(expected,model.electron.distribution.values), L"wrong distribution results", LINE_INFO());
		}

	};
}