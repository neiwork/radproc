#include <JetAGNResultChecks/JetAGNResultChecks.h>

#include <fparameters/parameters.h>
#include <JetAGN/ioutil.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using boost::property_tree::ptree;

int main()
{
	const std::string TEST_FOLDER = "../JetAGNUnitTests/";

	GlobalConfig = readConfig(TEST_FOLDER + "test-parameters.json");

	ptree results = JetAGNResultChecks::produceResults();

	boost::property_tree::write_json(TEST_FOLDER + JetAGNResultChecks::RESULTS_FILENAME, results);

	return 0;
}

