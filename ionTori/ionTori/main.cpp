#include <stdio.h>


#include "modelParameters.h"
#include "targetFields.h"
#include "State.h"
#include "write.h"
#include "radiativeLosses.h"

#include <inout/ioutil.h>


#include <fparticle/Particle.h>
#include <fparameters/parameters.h>
#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>
#include <fmath/physics.h>

#include <boost/property_tree/ptree.hpp>

int main()
{
	std::string folder{ prepareOutputfolder() };

	try {
		GlobalConfig = readConfig();
		prepareGlobalCfg();
		State model(GlobalConfig.get_child("model"));
		
		
		radiativeLosses(model, folder+"\\electronLosses.txt");
		
        
		//ParamSpaceValues psv(model.electron.ps);

		//psv.fill([&](const SpaceIterator& i){
		//	double E = i.val(DIM_E);
		//	double z = i.val(DIM_R);
		//	return frad(E, z);
		//});

		
		//injection(model.electron, model);
		
		//distribution(model.electron, model);

		//writeAllSpaceParam(getFileName(folder, "eDist"), model.electron.distribution);

		//processes(model, getFileName(folder, "luminosity"));

	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
//		throw;
	}

	return 0;
}




/*int main(int argc, char **argv)
{
	printf("hello world\n");
	return 0;
}*/
