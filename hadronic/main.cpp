#include "processes.h"
#include "injection.h"

#include "radiativeLosses.h"
#include "state.h"
#include "modelParameters.h"

#include <fparameters\parameters.h>
#include <inout\ioutil.h>

#include <boost\property_tree\ptree.hpp>

#include <iostream>

//#pragma comment(lib, "inout")

int main() {

	 

	std::string folder{ prepareOutputfolder() };

	try {
		GlobalConfig = readConfig();
		prepareGlobalCfg();
		State model(GlobalConfig.get_child("model"));

		radiativeLosses(model,folder+"\\protonLosses.txt");
		//writeEmax(folder + "\\Emax.txt", model.proton)

		injection(model.proton, model);
		distribution(model.proton, model);

		processes(model, folder + "\\luminosity.txt");

	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
		throw;
	}

	return 0;
}