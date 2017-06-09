#include "radiativeLosses.h"


#include "targetFields.h"
#include "write.h"

//#include "lossesAnisotropicIC.h"
#include <flosses\lossesSyn.h>
#include <flosses\nonThermalLosses.h>
#include <flosses\lossesIC.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>

void radiativeLosses(State& st, const std::string& filename, Vector& Gc)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	static const double starT = GlobalConfig.get<double>("starT");
	static const double starTIR = GlobalConfig.get<double>("IRstarT");
	
	double z_0 = st.electron.ps[DIM_R].first();

	

	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << "Log(E/eV)" //<< "\t" << "Log(R/pc)"
		<< "\t" << "Synchr"
		<< "\t" << "IC"
		//<< "\t" << "ICAux"
		<< "\t" << "Diff"
		<< "\t" << "Acc"
		<< "\t" << "Ad" << std::endl;
	
	double Emin = boltzmann*starT / 100.0;
	double EphminAux = boltzmann*starTIR / 100.0;
	
	double Reff = stagnationPoint(z_0);

	st.electron.ps.iterate([&](const SpaceIterator& i){
		const double B =  st.magf.get(i);
		double vel_lat = cLight*openingAngle;

		double E = i.val(DIM_E);
		double z = i.val(DIM_R);

		double gamma = Gc[i.coord[DIM_R]];

		double fmtE = log10(E / 1.6e-12);		

		double eSyn = lossesSyn(i.val(DIM_E), B, st.electron) / i.val(DIM_E);
		
		//VER como le paso el vector Gc a las perdidas
		double eIC2 = lossesIC(i.val(DIM_E), st.electron,
			[&E,&z,&gamma](double E){
			return starBlackBody(E,z,gamma);},   
				Emin, 1.0e4*Emin ) / i.val(DIM_E);
			
		double eIC_Aux = lossesIC(i.val(DIM_E), st.electron,
				[&E, &z, &gamma](double E) {
				return starIR(E, z, gamma); },
				EphminAux, 1.0e4*EphminAux) / i.val(DIM_E); 			

		double eDif  = diffusionRate(E, Reff, B);
		double eAcc = accelerationRate(E, B, accEfficiency);
		double eAdia = adiabaticLosses(E, z, vel_lat, gamma) / E;
		
		file << fmtE //<< "\t" << logR 
							<< "\t" << safeLog10(eSyn) 
							<< "\t" << safeLog10(eIC2)
							<< "\t" << safeLog10(eIC_Aux)
							<< "\t" << safeLog10(eDif)
							<< "\t" << safeLog10(eAcc)
							<< "\t" << safeLog10(eAdia) 
							<< std::endl;
	

	}); 


	file.close();
	



}
	//File protonLosses(files, out, std::string("protonLosses"));

	//st.proton.ps.iterate([&st, &files, &out](const SpaceIterator& i){

	//	double fmtE = log10(i.val(DIM_E) / 1.6e-12);

	//	double pSyn   = lossesSyn(i.val(DIM_E), st.proton) / i.val(DIM_E);
	//	double pp     = lossesHadronics(i.val(DIM_E), st.proton) / i.val(DIM_E);
	//	double pgamma = lossesPhotoHadronic(i.val(DIM_E), st.proton, st.tpf) / i.val(DIM_E);
	//	double pDif   = diffusionRate(i.val(DIM_E), st.proton);
	//	double pAcc   = accelerationRate(i.val(DIM_E), st.proton);

	//	out["protonLosses"]->file << fmtE	<< "\t" << safeLog10(pSyn)
	//										<< "\t" << safeLog10(pgamma)
	//										<< "\t" << safeLog10(pp)
	//										<< "\t" << safeLog10(pAcc)
	//										<< "\t" << safeLog10(pDif) << std::endl;


	//});

	