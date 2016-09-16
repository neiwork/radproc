#include "write.h"

#include "state.h"
#include "modelParameters.h"

#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>
#include <fparameters\parameters.h>
#include <boost/property_tree/ptree.hpp>

//namespace {
	double safeLog10( double x ) {
		return x>0.0 ? log10(x) : 0.0;
	}
//}

std::string dataName(std::string id) {
	return id + ".txt";
}

void generateViewScript(std::string filename) {
	std::ofstream file;
	file.open((filename+".bat").c_str(), std::ios::out);
	file << "@../plot-svg-and-view.bat " + filename.substr(filename.find("\\")+1);
	file.close();
}


void writeAllSpaceParam(const std::string& filename, const ParamSpaceValues& data)
{
	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	
	data.ps.iterate([&file, &data](const SpaceIterator& i){
		double logE = log10(i.val(DIM_E) / 1.6e-12);
		double lognu = log10(i.val(DIM_E) / planck);

		double logQ = log10(data.get(i)); //log10(salida.values(i));  // values(i));
//		salida.values(i);


		file << logE << '\t' << lognu << '\t' << 
			logQ << std::endl;
			//logQ << std::endl;
	});

	file.close();
	generateViewScript(filename);
}

double Llab(double Lint)
{
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	double boost = pow(Dlorentz, 4.0);
	return Lint*boost;
}

void dopplerBoost(const std::string& filename, const ParamSpaceValues& data)
{
	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	data.ps.iterate([&file, &data](const SpaceIterator& i) {

		double Elab = i.val(DIM_E)*Dlorentz; //Dlorentz=delta
		double logE = log10(Elab / 1.6e-12);
		double lognu = log10(Elab / planck);

		double logQ = log10(Llab(data.get(i)));


		file << logE << '\t' << lognu << '\t' <<
			logQ << std::endl;
		//logQ << std::endl;
	});

	file.close();
	generateViewScript(filename);
}


void readData(const std::string& archive, Matrix& data)
{
	std::ifstream fileReaded;
	fileReaded.open(archive.c_str(), std::ios::in);
	
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	double boost = pow(Dlorentz, 4.0);

	double logE, loglum, errLum;
	int j = 0;


	while (!fileReaded.eof())  //esto termina cuando llega al final
	{
		fileReaded >> logE >> loglum >> errLum;

		data[j][0] = pow(10.0, logE)*1.6e-12/Dlorentz; //paso del lab al ff
		data[j][1] = pow(10.0, loglum)/boost;
		data[j][2] = errLum / boost;

		j += 1;
	}

	fileReaded.close();
}


