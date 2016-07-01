#include "radiativeLosses.h"


#include "targetFields.h"
#include "write.h"

#include "lossesAnisotropicIC.h"
#include <flosses\lossesSyn.h>
#include <flosses\nonThermalLosses.h>
#include <flosses\lossesIC.h>
#include <fparameters\SpaceIterator.h>

#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>

void radiativeLosses(State& st)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle", 0.1);
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency", 0.1);
	
	std::vector<File*> files;
	OFM out;

	File electronLosses(files, out, std::string("electronLosses"));	

	out["electronLosses"]->file << "Log(E/eV)" << "\t" << "Log(R/pc)"
		<< "\t" << "Syn"
		<< "\t" << "IC"
		<< "\t" << "ICcmb"
		<< "\t" << "Dif"
		<< "\t" << "Acc"
		<< "\t" << "Adia" << std::endl;
	
	st.electron.ps.iterate([&st, &files, &out](const SpaceIterator& i){

		const double magf{ st.magf.get(i) };
		double fmtE = log10(i.val(DIM_E) / 1.6e-12);
		double logR = log10(i.val(DIM_R)/pc);
		//double logT = log10(i.val(DIM_T));

		double B = magf; // i.par.magneticField; VER por qué no funciona

		//double Reff = 10.0*stagnationPoint(i.val(DIM_R));
		double vel_lat = cLight*openingAngle;

		double E = i.val(DIM_E);

		double EphminS(0.0), EphminCMB(0.0);
		targetPhotonEnergies(EphminS, EphminCMB);

		double eSyn = lossesSyn(i.val(DIM_E), B, st.electron) / i.val(DIM_E);
		double eIC = lossesAnisotropicIC(i.val(DIM_E), st.electron, i.val(DIM_R)) / i.val(DIM_E);
		double eIC2 = lossesIC(i.val(DIM_E), st.electron,
			[E](double E){
			return cmbBlackBody(E);},   
				EphminCMB, 1.0e4*EphminCMB ) / i.val(DIM_E);
		
		   
		double eDif  = diffusionRate(i.val(DIM_E), i.val(DIM_R), B);
		double eAcc = accelerationRate(i.val(DIM_E), B, accEfficiency);
		double eAdia = adiabaticLosses(i.val(DIM_E), i.val(DIM_R), vel_lat) / i.val(DIM_E);
		
	out["electronLosses"]->file << fmtE << "\t" << logR 
											<< "\t" << safeLog10(eSyn) 
											<< "\t" << safeLog10(eIC)
											<< "\t" << safeLog10(eIC2)
											<< "\t" << safeLog10(eDif)
											<< "\t" << safeLog10(eAcc)
											<< "\t" << safeLog10(eAdia) 
											<< std::endl;
	

	}, { -1, -1, 0 }); //fijo t


	for (auto f : files) {
		f->file.close();
	}



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

	