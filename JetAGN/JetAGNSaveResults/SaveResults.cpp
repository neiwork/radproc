#include <JetAGNResultChecks/JetAGNResultChecks.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using boost::property_tree::ptree;

int main()
{
	ptree results = JetAGNResultChecks::produceResults();

	boost::property_tree::write_json("../JetAGNUnitTests/" + JetAGNResultChecks::RESULTS_FILENAME, results);

	return 0;
}

