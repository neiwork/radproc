#include "CppUnitTest.h"

#include <fparameters/Dimension.h>

#include <JetAGN/distribution.h>
#include <JetAGN/injection.h>
#include <JetAGN/State.h>
#include <JetAGN/modelParameters.h>
#include <JetAGN/checkPower.h>
#include <JetAGN/checks.h>

#include <boost/property_tree/ptree.hpp>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using boost::property_tree::ptree;

namespace JetAGNUnitTests
{
	/* Default parameter config for testing. */
	boost::property_tree::ptree config() {

		prepareGlobalCfg();

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
			State model(config());
			
			Assert::AreEqual<size_t>(4, model.electron.ps[0].size(), L"energy dimension size", LINE_INFO());
			Assert::AreEqual<size_t>(3, model.electron.ps[1].size(), L"energy dimension size", LINE_INFO());
			Assert::AreEqual<size_t>(3, model.electron.ps[2].size(), L"energy dimension size", LINE_INFO());
		}

		TEST_METHOD(TestEnergyDimensionValues11)
		{
			ptree modelcfg{ config() };
			modelcfg.put<int>("particle.default.dim.energy.samples", 11);
			State model(modelcfg);
			
			// check construction of standard 11-value energy vector
			std::vector<double> expected{ { 1.6000000000000001e-006, 1.6000000000000003e-005, 0.00016000000000000007, 0.0016000000000000009, 0.016000000000000011, 0.16000000000000014, 1.6000000000000016, 16.000000000000018, 160.00000000000020, 1600.0000000000023, 16000.000000000025 } };

			Assert::IsTrue(check_vec(expected, model.electron.ps[0].values), L"wrong energy dimension vector", LINE_INFO());
		}


		TEST_METHOD(TestInjection)
		{
			State model(config());

			injection(model.electron, model);

			double injectedPower = computeInjectedPower(model.electron.injection, 0);
			Assert::AreEqual(1.1372460991317841e+044, injectedPower, L"wrong injected power", LINE_INFO());

		}

		TEST_METHOD(TestDistribution)
		{
			State model(config());

			injection(model.electron, model);
			distribution(model.electron, model);

			// HARDCODED EXPECTED RESULTS IN THE PSV (36 values)
			std::vector<double> expected{ { 0.0030372331498481585, 6.7271865939765752e-008, 2.4502560741810692e-019, 0.00000000000000000, 2.5851187560565028, 7.1309890085422364e-005, 3.9832401155735729e-014, 0.00000000000000000, 4.4527859555133960, 0.00016157876475405818, 1.0821169757914407e-009, 0.00000000000000000, 7.8164621086299211e-008, 2.0661402204719002e-011, 2.4439381233911514e-017, 0.00000000000000000, 0.0087534231598507487, 5.0787750129479747e-006, 1.1085089364559843e-010, 0.00000000000000000, 209.56724283996866, 0.12398380704075215, 3.3274439247937427e-005, 0.00000000000000000, 7.7963074711314374e-006, 2.0608127070027726e-009, 2.4376364633483944e-015, 0.00000000000000000, 24.360181344700159, 0.014133887554961632, 3.0849058010299499e-007, 0.00000000000000000, 6444057.5919388374, 3812.4222122133679, 1.0231686831921796, 0.00000000000000000 } };

			Assert::IsTrue(check_vec(expected,model.electron.distribution.values), L"wrong distribution results", LINE_INFO());
		}

	};
}