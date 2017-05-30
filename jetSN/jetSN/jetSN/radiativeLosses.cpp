#include "radiativeLosses.h"


#include "targetFields.h"
#include "write.h"

//#include "lossesAnisotropicIC.h"
#include <flosses\lossesSyn.h>
#include <flosses\nonThermalLosses.h>
#include <flosses\lossesIC.h>
#include <fparameters\SpaceIterator.h>

#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>

void radiativeLosses(State& st, const std::string& filename)
{
	//static const double Bfield = GlobalConfig.get<double>("Bfield");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	static const double starT = GlobalConfig.get<double>("starT");
	
	double r = GlobalConfig.get<double>("z_int")*pc;
	
	//std::vector<File*> files;
	//OFM out;

	//File electronLosses(files, out, std::string("electronLosses"));	

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
		
	double B = computeMagField(r);
	double vel_lat = cLight*openingAngle;
	double Reff = stagnationPoint(r);

	st.electron.ps.iterate([&](const SpaceIterator& i){
		//const double magf{ st.magf.get(i) };

		double E = i.val(DIM_E);
		double fmtE = log10(E / 1.6e-12);		

		double eSyn = lossesSyn(i.val(DIM_E), B, st.electron) / i.val(DIM_E);
		
		double eIC2 = lossesIC(i.val(DIM_E), st.electron,
			[&E,&r](double E){
			return starBlackBody(E,r);},   
				Emin, 1.0e4*Emin ) / i.val(DIM_E);
			
		/*double eIC_Aux = lossesIC(i.val(DIM_E), st.electron,
				[&E, &r](double E) {
				return starIR(E, r); },
				EphminAux, 1.0e4*EphminAux) / i.val(DIM_E); */
				

		double eDif  = diffusionRate(E, Reff, B);
		double eAcc = accelerationRate(E, B, accEfficiency);
		double eAdia = adiabaticLosses(E, r, vel_lat) / E;
		
		file << fmtE //<< "\t" << logR 
							<< "\t" << safeLog10(eSyn) 
							<< "\t" << safeLog10(eIC2)
							//<< "\t" << safeLog10(eIC_Aux)
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

	