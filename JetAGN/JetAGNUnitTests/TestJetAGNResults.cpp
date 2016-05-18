#include "CppUnitTest.h"

#include <JetAGN/distribution.h>
#include <JetAGN/injection.h>
#include <JetAGN/State.h>
#include <JetAGN/modelParameters.h>
#include <JetAGN/checkPower.h>
#include <JetAGN/checks.h>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace JetAGNUnitTests
{

	TEST_CLASS(JetAGNResults)
	{
	public:

		TEST_METHOD(Injection)
		{
			setParameters();

			const double energyPointsForTesting = 4;
			ParticleCfg<Electron>::config.nE = energyPointsForTesting;
			ParticleCfg<Photon>::config.nE = energyPointsForTesting;

			State model(3);

			injection(model.electron, model);

			double injectedPower = computeInjectedPower(model.electron.injection, 0);
			Assert::AreEqual(1.8828708877129821e+044, injectedPower, L"wrong injected power", LINE_INFO());

		}

		TEST_METHOD(Distribution)
		{
			setParameters();

			// OVERRIDE DISCRETIZATION TO MAKE SURE TEST RUNS FAST
			// >>
			const double energyPointsForTesting = 3;
			ParticleCfg<Electron>::config.nE = energyPointsForTesting;
			ParticleCfg<Photon>::config.nE = energyPointsForTesting;
			State model(2);
			// <<

			injection(model.electron, model);
			distribution(model.electron, model);

			// HARDCODED EXPECTED RESULTS IN THE PSV (36 values)
			std::vector<double> expected{ { 0.0030372331498481585, 6.7271865939765752e-008, 2.4502560741810692e-019, 0.00000000000000000, 2.5851187560565028, 7.1309890085422364e-005, 3.9832401155735729e-014, 0.00000000000000000, 4.4527859555133960, 0.00016157876475405818, 1.0821169757914407e-009, 0.00000000000000000, 7.8164621086299211e-008, 2.0661402204719002e-011, 2.4439381233911514e-017, 0.00000000000000000, 0.0087534231598507487, 5.0787750129479747e-006, 1.1085089364559843e-010, 0.00000000000000000, 209.56724283996866, 0.12398380704075215, 3.3274439247937427e-005, 0.00000000000000000, 7.7963074711314374e-006, 2.0608127070027726e-009, 2.4376364633483944e-015, 0.00000000000000000, 24.360181344700159, 0.014133887554961632, 3.0849058010299499e-007, 0.00000000000000000, 6444057.5919388374, 3812.4222122133679, 1.0231686831921796, 0.00000000000000000 } };

			Assert::IsTrue(check_vec(expected,model.electron.distribution.values), L"wrong distribution results", LINE_INFO());
		}

	};
}