#include <JetAGN/checks.h>

#include <JetAGN/distribution.h>
#include <JetAGN/injection.h>
#include <JetAGN/checkPower.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>
#include <JetAGN/State.h>
#include <JetAGN/modelParameters.h>
#include <JetAGN/checks.h>

#include <boost/property_tree/ptree.hpp>

namespace JetAGNResultChecks {

	using boost::property_tree::ptree;

	/* Default parameter config for testing. */
	boost::property_tree::ptree config() {
		prepareGlobalCfg();
		return GlobalConfig.get_child("model");
	}

	std::vector<double> computeDistribution() {
		State model(config());
		injection(model.electron, model);
		distribution(model.electron, model);
		return std::move(model.electron.distribution.values);
	}

	std::vector<double> computeInjection() {
		State model(config());
		injection(model.electron, model);
		double injectedPower = computeInjectedPower(model.electron.injection, 0);
		return std::move(std::vector<double>({ injectedPower }));
	}

	std::vector<double> computeEnergyDimValues() {
		State model(config());
		return std::move(model.electron.ps[0].values);
	}

	ptree vec2ptree(std::vector<double> v) {
		ptree output;
		for (double n : v) {
			ptree child;
			child.put("", n);
			output.push_back(std::make_pair("", std::move(child)));
		}
		return std::move(output);
	}

	ptree __declspec(dllexport) produceResults() {
		ptree results;
		results.add_child("distribution", vec2ptree(computeDistribution()));
		results.add_child("injection", vec2ptree(computeInjection()));
		results.add_child("energy-dim-values", vec2ptree(computeEnergyDimValues()));
		return std::move(results);
	};

}