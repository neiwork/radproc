#include <boost/property_tree/ptree.hpp>
#include <vector>

namespace JetAGNResultChecks {
	const std::string RESULTS_FILENAME = "jet-agn-test-results.json";
	using boost::property_tree::ptree;
	ptree config();

	std::vector<double> computeDistribution();

	std::vector<double> computeInjection();

	std::vector<double> computeEnergyDimValues();

	boost::property_tree::ptree produceResults();

	template <typename T>
	std::vector<T> as_vector(ptree const& pt, ptree::key_type const& key)
	{
		std::vector<T> r;
		for (auto& item : pt.get_child(key))
			r.push_back(item.second.get_value<T>());
		return r;
	}
}