#include "CppUnitTest.h"


#include <JetAGNResultChecks/JetAGNResultChecks.h>
#include <JetAGN/modelParameters.h>
#include <JetAGN/checks.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>
#include <JetAGN/ioutil.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using boost::property_tree::ptree;
using JetAGNResultChecks::as_vector;

namespace JetAGNUnitTests
{

#define STRINGIFY(x) #x
#define EXPAND(x) STRINGIFY(x)

	static const std::string project_dir_full = EXPAND(UNITTESTPRJ);

	static ptree EXPECTED;

	TEST_MODULE_INITIALIZE(loadConfig) {
		std::string project_dir = project_dir_full;
		project_dir.erase(0, 1); // erase the first quote
		project_dir.erase(project_dir_full.size() - 4); // erase the last quote and the dot
		boost::property_tree::read_json(project_dir + "\\" + JetAGNResultChecks::RESULTS_FILENAME, EXPECTED);
		GlobalConfig = readConfig(project_dir + "\\test-parameters.json");
	}

	TEST_CLASS(JetAGNResults)
	{

	public:

		TEST_METHOD(TestStateConfiguration)
		{
			State model(JetAGNResultChecks::config());
			
			Assert::AreEqual<int>(4, model.electron.ps[0].size(), L"energy dimension size", LINE_INFO());
			Assert::AreEqual<int>(3, model.electron.ps[1].size(), L"energy dimension size", LINE_INFO());
			Assert::AreEqual<int>(3, model.electron.ps[2].size(), L"energy dimension size", LINE_INFO());
		}

		TEST_METHOD(TestEnergyDimensionValues)
		{
			std::vector<double> results = JetAGNResultChecks::computeEnergyDimValues();
			std::vector<double> expected = as_vector<double>(EXPECTED, "energy-dim-values");
			Assert::IsTrue(check_vec(expected, results), L"energy dimension vector", LINE_INFO());
		}


		TEST_METHOD(TestInjection)
		{
			std::vector<double> results = JetAGNResultChecks::computeInjection();
			std::vector<double> expected = as_vector<double>(EXPECTED, "injection");
			Assert::IsTrue(check_vec(expected, results), L"injected power", LINE_INFO());

		}

		TEST_METHOD(TestDistribution)
		{
			std::vector<double> results = JetAGNResultChecks::computeDistribution();
			std::vector<double> expected = as_vector<double>(EXPECTED, "distribution");
			Assert::IsTrue(check_vec(expected,results), L"distribution results", LINE_INFO());
		}

	};

}