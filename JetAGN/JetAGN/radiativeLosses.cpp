#include "radiativeLosses.h"


#include "targetFields.h"
#include "write.h"

#include "lossesAnisotropicIC.h"
#include <flosses\lossesSyn.h>
#include <flosses\nonThermalLosses.h>
#include <flosses\lossesIC.h>
#include <fparameters\SpaceIterator.h>

//#include "losses.h"
//#include <flosses\lossesHadronics.h>
//#include <flosses\lossesPhotoHadronic.h>
//#include <flosses\lossesBrem.h>
//#include <fmath\interpolation.h>




void radiativeLosses(State& st)
{
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


		double fmtE = log10(i.par.E / 1.6e-12);
		double logR = log10(i.par.R/pc);
		//double logT = log10(i.par.T);

		double B = parameters.magneticField; // i.par.magneticField; VER por qué no funciona

		double Reff = 10.0*stagnationPoint(i.par.R);
		double vel_lat = cLight*parameters.openingAngle;

		double E = i.par.E;

		double EphminS(0.0), EphminCMB(0.0);
		targetPhotonEnergies(EphminS, EphminCMB);

		double eSyn = lossesSyn(i.par.E, B, st.electron) / i.par.E;
		double eIC = lossesAnisotropicIC(i.par.E, st.electron, i.par.R) / i.par.E;
		double eIC2 = lossesIC(i.par.E, st.electron,
			[E](double E){
			return cmbBlackBody(E);},   
				EphminCMB, 1.0e4*EphminCMB ) / i.par.E;
		
		   
		double eDif  = diffusionRate(i.par.E, i.par.R, B);
		double eAcc = accelerationRate(i.par.E, B, parameters.accEfficiency);
		double eAdia = adiabaticLosses(i.par.E, i.par.R, vel_lat) / i.par.E;
		
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

	//	double fmtE = log10(i.par.E / 1.6e-12);

	//	double pSyn   = lossesSyn(i.par.E, st.proton) / i.par.E;
	//	double pp     = lossesHadronics(i.par.E, st.proton) / i.par.E;
	//	double pgamma = lossesPhotoHadronic(i.par.E, st.proton, st.tpf) / i.par.E;
	//	double pDif   = diffusionRate(i.par.E, st.proton);
	//	double pAcc   = accelerationRate(i.par.E, st.proton);

	//	out["protonLosses"]->file << fmtE	<< "\t" << safeLog10(pSyn)
	//										<< "\t" << safeLog10(pgamma)
	//										<< "\t" << safeLog10(pp)
	//										<< "\t" << safeLog10(pAcc)
	//										<< "\t" << safeLog10(pDif) << std::endl;


	//});

	