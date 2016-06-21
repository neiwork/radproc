#include "maximumEnergy.h"


#include "modelParameters.h"

//#include "state.h"

#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>
//#include <fparameters\parameters.h>
//#include <boost/property_tree/ptree.hpp>


//std::string dataName(std::string id) {
//	return id + ".txt";
//}


void writeEmax(const std::string& filename, Particle& p)
{

	std::ofstream file;
	//file.open(dataName(filename).c_str(), std::ios::out);
	file.open((filename + ".txt").c_str(), std::ios::out);

	// version acotada
	//double time = log10(data.ps[1][t]);

	//for (int r_ix = 0; r_ix < p.ps[1].size(); r_ix++) {

	p.ps.iterate([&](const SpaceIterator& i){

		double z = i.val(DIM_R);
		double logR = ( z / pc);

		double B = computeMagField(z);

		double Emax = eEmax(z, B);
		double logE = log10(Emax/1.6e-12);

		file << logR << '\t' << logE << std::endl;
		}, { 0, -1, 0 });

	file.close();
	//generateViewScript(filename);
}
