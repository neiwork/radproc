#include "write.h"

#include "state.h"
#include "modelParameters.h"
#include <fparameters\parameters.h>

//namespace {
	double safeLog10( double x ) {
		return x>0.0 ? log10(x) : 0.0;
	}
//}




void writeAllSpaceParam(const std::string& filename, const ParamSpaceValues& data)
{
	
	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);
	//const ParamSpace* a = &(data.ps);
	
	// a.iterate;  no me deja hacer esta operacion
	
	data.ps.iterate([&file, &data](const SpaceIterator& i){

		double logE = log10(i.par.E / 1.6e-12);
		double logR = log10(i.par.R);
		double logT = log10(i.par.T);
		double logQ = log10(data.get(i)); //log10(salida.values(i));  // values(i));
//		salida.values(i);


		file << logE << '\t' << logR << '\t' << logT << '\t' << 
			logQ << std::endl;
			//logQ << std::endl;
	});

	file.close();
}

void writeEandTParamSpace(const std::string& filename, const ParamSpaceValues& data, int r)
{

	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	// version acotada
	double logR = log10(data.ps[1][r]);
	

	file << "log(r)=" << logR << '\t' << std::endl;
	data.ps.iterate([&file, &data](const SpaceIterator& i){

		double logE = log10(i.par.E / 1.6e-12);
		double time = i.par.T;
		double logQ = log10(data.get(i));

		file << logE << '\t' << time << '\t' << logQ << std::endl;
		;
	}, { -1, r, -1 });  //el -1 indica que las E se recorren, no quedan fijas
	//las otras dos dimensiones quedan fijas en las posiciones r y t (recordar que la primera es 0 )

	file.close();
}

void writeEnergyFunction(const std::string& filename, const ParamSpaceValues& data, int r, int t)
{

	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	// version acotada
	double logR = log10(data.ps[1][r]);
	// version larga
	double logT = log10(data.ps.dimensions[2]->values[t]);

	file << "log(r)=" << logR << '\t' << "log(t)=" << logT << std::endl;
	data.ps.iterate([&file, &data](const SpaceIterator& i){

		double logE = log10(i.par.E / 1.6e-12);
		double logQ = log10(data.get(i));

		file << logE << '\t' << logQ << std::endl;
		;
	}, { -1, r, t });  //el -1 indica que las E se recorren, no quedan fijas
	//las otras dos dimensiones quedan fijas en las posiciones r y t (recordar que la primera es 0 )

	file.close();
}


void write(const std::string& archive, Vector salida, Particle particle)
{
	std::ofstream fileWritten;
	fileWritten.open(archive.c_str(), std::ios::out);

	for (size_t i = 0; i < particle.eDim()->size(); ++i){

		double E_eV = (particle.eDim()->values[i] / 1.6e-12);

		fileWritten << log10(E_eV);

		double result = salida[i];

		fileWritten << "\t" << log10(result) << std::endl;

	}

	fileWritten.close();
}



