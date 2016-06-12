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
			Assert::AreEqual(7.5643007696355601e+033, injectedPower, L"wrong injected power", LINE_INFO());

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
			std::vector<double> expected{ { 1.6039693002058387e-007, 3.7311183244389524e-012, 1.4369049899920965e-023, 0.00000000000000000, 5.0721962856368329e-009, 1.1798832124818175e-013, 4.5438925496364766e-025, 0.00000000000000000, 1.6039693002058388e-010, 3.7311183244389527e-015, 1.4369049899920967e-026, 0.00000000000000000, 1.6040151383267041e-007, 3.7323299721618590e-012, 1.4475690006056270e-021, 0.00000000000000000, 5.0724410933991829e-009, 1.1808457317714854e-013, 1.3103108807146532e-021, 0.00000000000000000, 4.7597059198617138e-008, 1.5070319044837016e-012, 4.8327517200303899e-022, 0.00000000000000000, 1.6085871311183251e-007, 3.8531823227930523e-012, 1.4439801548170083e-019, 0.00000000000000000, 5.3747482066771741e-009, 2.8891767214281165e-013, 3.6510733020081490e-018, 0.00000000000000000, 2.6355740401765808e-010, 6.4507673785987928e-014, 1.4975873293337995e-017, 0.00000000000000000 } };

			Assert::IsTrue(check_vec(expected,model.electron.distribution.values), L"wrong distribution results", LINE_INFO());
		}

	};
}