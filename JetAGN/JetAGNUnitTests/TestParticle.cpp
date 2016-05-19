#include "CppUnitTest.h"

#include <fparticle/Particle.h>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace JetAGNUnitTests
{
	TEST_CLASS(TestParticle)
	{
	public:
		
		TEST_METHOD(TestParticleConfig)
		{
			boost::property_tree::ptree cfg;
			cfg.add<double>("mass", 1);
			cfg.add<double>("dim.energy.min", 2);
			cfg.add<double>("dim.energy.max", 3);

			Particle p("myparticle");
			p.configure(cfg);

			Assert::AreEqual(double{ 1.0 }, p.mass, L"mass not set", LINE_INFO());
			Assert::AreEqual(double{ 2.0 }, p.logEmin, L"emin not set", LINE_INFO());
			Assert::AreEqual(double{ 3.0 }, p.logEmax, L"emax not set", LINE_INFO());

		}

	};
}