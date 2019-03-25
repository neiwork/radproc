#include "write.h"

#include "state.h"
#include "state.h"
#include "modelParameters.h"


#include "targetFieldM87.h"
//#include "lossesAnisotropicIC.h"
#include <flosses\lossesSyn.h>
#include <flosses\nonThermalLosses.h>


#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>
#include <fparameters\parameters.h>
#include <boost/property_tree/ptree.hpp>


std::string dataName(std::string id) {
	return id + ".txt";
}

//void generateViewScript(std::string path) {
//	std::string filename = path.substr(path.find("\\") + 1);
//	std::string folder = path.substr(0,path.find("\\"));
//	std::ofstream file;
//	file.open((folder+"/plots/plot-"+filename+".bat").c_str(), std::ios::out);
//	file << "@../../plot-svg-and-view.bat " + filename;
//	file.close();
//}


void writeAllSpaceParam(const std::string& filename, const ParamSpaceValues& data, 
	const Vector& Gc, const Vector& tobs)
{
	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);
	
	file << "log(E / eV)" << '\t' <<
		"log(z / pc)" << '\t' <<
		"t / yr" << '\t' <<
		// logT << '\t' << 
		" logQ"  << std::endl;
	//logQ << std::endl;
	
	const double z_0 = data.ps[DIM_R].first();
	double t0 = z_0 / cLight;

	data.ps.iterate([&](const SpaceIterator& i){
		double logE = log10(i.val(DIM_E) / 1.6e-12);
		double z = i.val(DIM_R);
		double logR = log10(z/pc);
		double logQ = log10(data.get(i));

		double gamma = Gc[i.coord[DIM_R]];
		
		double tlab = z / cLight - t0;
		double T = tobs[i.coord[DIM_R]];

		file << logE << '\t' << 
			logR << '\t' <<
			 T/yr << '\t' << 
			logQ << std::endl;
			//logQ << std::endl;
	});

	file.close();
	//generateViewScript(filename);
}

void writeEnergyFunction(const std::string& filename, const ParamSpaceValues& data, int r)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	// version acotada
	//double logR = log10(data.ps[1][r]);

	data.ps.iterate([&file, &data](const SpaceIterator& i) {

		double logE = log10(i.val(DIM_E) / 1.6e-12);
		double logQ = safeLog10(data.get(i));

		file << logE << '\t' << logQ << std::endl;
		;
	}, { -1, r});  
					  

	file.close();
	
}




/*

void writeEvol(const std::string& archive, const std::string& filename, 
	State& st)
{

	static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double IRstarT = GlobalConfig.get<double>("IRstarT");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");


	std::ifstream fileReaded;
	fileReaded.open(archive.c_str(), std::ios::in);

	if (!fileReaded.is_open())
	{
		std::cout << "el archivo para leer no se puedo abrir";
	}
	
	//if (fileReaded.good()) {
	//	std::cout << "el archivo esta good";
	//}

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	int j = 0;

	double tobs, t, z_norm, g, ratio, rho, Lnt;
	
	
	while (!fileReaded.eof())  //esto termina cuando llega al final
		{
			fileReaded >> tobs;
			fileReaded >> t;
			fileReaded >> z_norm;
			fileReaded >> g;
			fileReaded >> ratio;
			fileReaded >> rho;
			fileReaded >> Lnt;

			if (j==0 || j % 100 == 0) {
				double z = z_norm*pc;
				double Gc = g*Gj;
				double Rc = ratio*z*openingAngle;

				double E = P2(electronMass*cLight2) / (2.7*boltzmann*IRstarT) / Gc;

				double vel_lat = cLight*openingAngle;

				double tad = adiabaticLosses(E, z, vel_lat, Gc) / E; //FF

				double B = computeMagField(z, Gc);

				double eIC_Aux = lossesAnisotropicIC(E, st.electron, z, Gc,
					[&](double E, double z, double r) {
					return nph_ICani(E, z, r, Gc, "IR"); },
					IRstarT) / E;

				double eSyn = lossesSyn(E, B, st.electron) / E;

				double frad = eIC_Aux / (tad + eIC_Aux + eSyn);

				file << z_norm
					<< '\t' << frad
					<< '\t' << eIC_Aux
					<< '\t' << eSyn
					<< '\t' << B << std::endl;
			}
			j += 1;
		}

	fileReaded.close();
	
	file.close();

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




void writeEnt(const std::string& filename, const ParamSpaceValues& data)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gamma = GlobalConfig.get<double>("Gamma");

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

//<<<<<<< HEAD
	file << "log(z/pc)" << '\t';
//=======
//	for (int i = 0; i < particle.eDim()->size(); ++i){
//>>>>>>> globals

	for (size_t t_ix = 0; t_ix < data.ps[2].size(); t_ix++) {
		double time = data.ps[2][t_ix];
		file << "t=" << log10(time) << '\t';
	}

	const double Emin = data.ps[DIM_E].first();

	const double RMIN = data.ps[DIM_R].first();
	const double RMAX = data.ps[DIM_R].last();
	const int N_R = data.ps[DIM_R].size() - 1;

	double z_int = pow((RMAX / RMIN), (1.0 / N_R));

	for (size_t r_ix = 0; r_ix < data.ps[1].size()-1; r_ix++) {

		file << std::endl;

		double z = data.ps[1][r_ix];
		double logR =  (z / pc);

		double dz = z * (z_int - 1);

		//volumen de la celda i
		double vol_i = pi*P2(jetRadius(z, openingAngle))*dz;;

		// [av] ver: no se usaba?
		//double Emax = eEmax(z, magneticField);

		file << logR << '\t';

		for (size_t t_ix = 0; t_ix < data.ps[2].size(); t_ix++) {

			double sum = 0.0;

			for (size_t E_ix = 0; E_ix < data.ps[0].size()-1; E_ix++) {

				double dE = data.ps[0][E_ix + 1] - data.ps[0][E_ix];
				double E = data.ps[0][E_ix];
				double T = data.ps[2][t_ix];

				double dist = data.interpolate({ { DIM_E, E }, { DIM_R, z }, { DIM_T, T } });

				sum = sum + dist*E*dE;// *vol_i;
			}

			file << log10(sum*Gamma) << '\t';
			       //multiplico nuevamente por Gamma para obtener la Ent en el sist lab
		}

	}
	file.close();
	generateViewScript(filename);
}
*/


