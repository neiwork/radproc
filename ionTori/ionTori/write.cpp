#include "write.h"

#include "State.h"
#include "modelParameters.h"

#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>

//namespace {
	double safeLog10( double x ) {
		return x>0.0 ? log10(x) : 0.0;
	}
//}

std::string dataName(std::string id) {
	return id + ".txt";
}

void generateViewScript(std::string path) {
	std::string filename = path.substr(path.find("\\") + 1);
	std::string folder = path.substr(0,path.find("\\"));
	std::ofstream file;
	file.open((folder+"/plots/plot-"+filename+".bat").c_str(), std::ios::out);
	file << "@../../plot-svg-and-view.bat " + filename;
	file.close();
}


void writeAllSpaceParam(const std::string& filename, const ParamSpaceValues& data)
{
	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);
	//const ParamSpace* a = &(data.ps);
	
	// a.iterate;  no me deja hacer esta operacion
	
	data.ps.iterate([&file, &data](const SpaceIterator& i){
		double logE = log10(i.val(DIM_E) / 1.6e-12);
		double logR = log10(i.val(DIM_R));
		double logT = log10(i.val(DIM_THETA));
		double logQ = log10(data.get(i)); //log10(salida.values(i));  // values(i));
//		salida.values(i);


		file << logE << '\t' << logR << '\t' << logT << '\t' << 
			logQ << std::endl;
			//logQ << std::endl;
	});

	file.close();
	generateViewScript(filename);
}

void writeEandRParamSpace(const std::string& filename, const ParamSpaceValues& data, int t)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	// version acotada
	//double time = log10(data.ps[1][t]);

	for (int r_ix = 0; r_ix < data.ps[1].size(); r_ix++) {

		double logR = log10(data.ps[1][r_ix] / pc);
		
		data.ps.iterate([&](const SpaceIterator& i){

			double logE = log10(i.val(DIM_E) / 1.6e-12);

			double logQ = safeLog10(data.get(i));

			file << logE << '\t' << logR << '\t' << logQ << std::endl;;
			
		}, { -1, r_ix, t });  
	}
	file.close();
	generateViewScript(filename);
}

void writeEandTParamSpace(const std::string& filename, const ParamSpaceValues& data, int r)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	// version acotada
	double logR = log10(data.ps[1][r]);
	
	file << "log(r)=" << logR << '\t' ;

	for (size_t t_ix = 0; t_ix < data.ps[2].size(); t_ix++) {
		double time = data.ps[2][t_ix];
		file << "t=" << log10(time) << '\t';
	}


	for (int E_ix = 0; E_ix < data.ps[0].size(); E_ix++) {

		file << std::endl;

		double logE = log10(data.ps[0][E_ix]/1.6e-12);

		file << logE << '\t';

		data.ps.iterate([&file, &data](const SpaceIterator& i){

			//double logR = log10(i.val(DIM_R));
			//double time = i.val(DIM_T);
			double logQ = safeLog10(data.get(i));

			file << logQ << '\t';
			;
		}, { E_ix, r, -1 });  //el -1 indica que las E se recorren, no quedan fijas
		//las otras dos dimensiones quedan fijas en las posiciones r y t (recordar que la primera es 0 )
	}
	file.close();
	generateViewScript(filename);
}


void writeRandTParamSpace(const std::string& filename, const ParamSpaceValues& data, int E)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	// version acotada
	double logE = log10(data.ps[0][E]/1.6e-12);


	file << "log(E)=" << logE << '\t' ;

	for (size_t t_ix = 0; t_ix < data.ps[2].size(); t_ix++) {
		double time = data.ps[2][t_ix];
		file << "t=" << log10(time) << '\t';
	}


	for (int r_ix = 0; r_ix < data.ps[1].size(); r_ix++) {

		file << std::endl;

		double logR = data.ps[1][r_ix] / pc;

		file << logR << '\t';

		data.ps.iterate([&file, &data](const SpaceIterator& i){

			double logQ = safeLog10(data.get(i));

			file << logQ << '\t';
			;
		}, { E, r_ix, -1 });  //el -1 indica que las E se recorren, no quedan fijas
		//las otras dos dimensiones quedan fijas en las posiciones r y t (recordar que la primera es 0 )
	}
	file.close();
	generateViewScript(filename);
}

void writeEnergyFunction(const std::string& filename, const ParamSpaceValues& data, int r, int t)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	// version acotada
	double logR = log10(data.ps[1][r]);
	// version larga
	double logT = log10(data.ps.dimensions[2]->values[t]);

	file << "log(r)=" << logR << '\t' << "log(t)=" << logT << std::endl;
	data.ps.iterate([&file, &data](const SpaceIterator& i){

		double logE = log10(i.val(DIM_E) / 1.6e-12);
		double logQ = safeLog10(data.get(i));

		file << logE << '\t' << logQ << std::endl;
		;
	}, { -1, r, t });  //el -1 indica que las E se recorren, no quedan fijas
	//las otras dos dimensiones quedan fijas en las posiciones r y t (recordar que la primera es 0 )

	file.close();
	generateViewScript(filename);
}

