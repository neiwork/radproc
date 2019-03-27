#include "radiativeLosses.h"


#include "targetFieldM87.h"


//#include "lossesAnisotropicIC.h"
#include <flosses\lossesSyn.h>
#include <flosses\lossesPhotoHadronic.h>
#include <flosses\nonThermalLosses.h>
#include <flosses\lossesIC.h>
#include <fparameters\SpaceIterator.h>

#include <fparameters\parameters.h>

#include <boost\property_tree\ptree.hpp>

void radiativeLosses(State& st, const std::string& filename)
{
	static const std::string id = GlobalConfig.get<std::string>("id");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	static const double starT = GlobalConfig.get<double>("starT");
	static const double Gamma = GlobalConfig.get<double>("Gb");
	
	double z = GlobalConfig.get<double>("z_peak")*pc;

	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << "Log(E/eV)" //<< "\t" << "Log(R/pc)"
		<< "\t" << "Synchr"
		<< "\t" << "pg"
		<< "\t" << "IC"
		<< "\t" << "IC-cmb"
		<< "\t" << "Acc"
		<< "\t" << "Ad" << std::endl;
	
	double Emin = boltzmann*starT / 100.0;
	double EphminAux;


	static const double cmbT = GlobalConfig.get<double>("cmbT");
	EphminAux = boltzmann*cmbT / 100.0;


	double B = computeMagField(z);

	double Emax = eEmax(z, B);

	st.proton.ps.iterate([&](const SpaceIterator& i){

		//const double magf{ st.magf.get(i) };
		double fmtE = log10(i.val(DIM_E) / 1.6e-12);

		//double r = i.val(DIM_R);
		//double logR = log10(r/pc);

		//double B = magf; 

		double vel_lat = cLight*openingAngle;

		double E = i.val(DIM_E);

		double eSyn = lossesSyn(i.val(DIM_E), B, st.proton) / i.val(DIM_E);


		double pGam = lossesPhotoHadronic(i.val(DIM_E), st.proton, [&E, &z](double E) {
			return starBlackBody(E, z); },
			Emin, 1.0e4*Emin) / i.val(DIM_E);

		
		double eIC2 = lossesIC(i.val(DIM_E), st.proton,
			[&E,&z](double E){
			return starBlackBody(E,z);},   
				Emin, 1.0e4*Emin ) / i.val(DIM_E);
			
		double eIC_Aux = lossesIC(i.val(DIM_E), st.proton,
				[&E, &z](double E) {
				return cmbBlackBody(E, z); },
				EphminAux, 1.0e4*EphminAux) / i.val(DIM_E);

				
		//double Reff = 10.0*stagnationPoint(i.val(DIM_R));
		//double eDif  = diffusionRate(i.val(DIM_E), Reff, B);
		double eAcc = accelerationRate(i.val(DIM_E), B, accEfficiency);
		double eAdia = adiabaticLosses(i.val(DIM_E), z, vel_lat, Gamma) / i.val(DIM_E);
		
		file << fmtE //<< "\t" << logR 
							<< "\t" << safeLog10(eSyn)
							<< "\t" << safeLog10(pGam)
							<< "\t" << safeLog10(eIC2)
							<< "\t" << safeLog10(eIC_Aux)
							//<< "\t" << safeLog10(eDif)
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

	