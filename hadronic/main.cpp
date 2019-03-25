#include "radiativeLosses.h"
#include "state.h"
#include "modelParameters.h"

#include <fparameters\parameters.h>
#include <inout\ioutil.h>

#include <boost\property_tree\ptree.hpp>

#include <iostream>

int main() {



	std::string folder{ prepareOutputfolder() };

	try {
		GlobalConfig = readConfig();
		prepareGlobalCfg();
		State model(GlobalConfig.get_child("model"));

		radiativeLosses(model,folder+"\\protonLosses.txt");

	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
		throw;
	}

	return 0;
}