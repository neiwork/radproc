#include "radiativeLosses.h"



#include "write.h"

#include <flosses\lossesSyn.h>
#include <flosses\nonThermalLosses.h>

//#include "losses.h"
//#include <flosses\lossesIC.h>
//#include <flosses\lossesHadronics.h>
//#include <flosses\lossesPhotoHadronic.h>
//#include <flosses\lossesBrem.h>
//#include <fmath\interpolation.h>




void radiativeLosses(State& st)
{
	std::vector<File*> files;
	OFM out;

	File electronLosses(files, out, std::string("electronLosses"));	

	//st.electron.ps.iterate([&st, &files, &out](const SpaceIterator& i){

		st.electron.ps.iterate([&st, &files, &out](const SpaceIterator& i){


		double fmtE = log10(i.par.E / 1.6e-12);
		double logR = log10(i.par.R);
		double logT = log10(i.par.T);

		double Reff = 10.0*(i.par.R*openingAngle); //revisar porque esto deberia ser el stagnation point


		double eSyn  = lossesSyn(i.par.E, st.electron) / i.par.E;		
		double eDif = diffusionRate(i.par.E, i.par.R);//, magneticField); // (double E, double radius, double));
		double eAcc = accelerationRate(i.par.E);//, magneticField, accEfficiency);
		double eAdia = adiabaticLosses(i.par.E, i.par.R, cLight) / i.par.E;

		//debug_dump();


		//double eIC   = lossesIC(i.par.E, st.electron, st.tpf) / i.par.E;
		//double eBrem = lossesBremss(i.par.E, st.electron) / i.par.E;
		
		out["electronLosses"]->file << fmtE << "\t" << logR << "\t" << logT << "\t" << safeLog10(eSyn)
										//	<< "\t" << safeLog10(eIC)
										//	<< "\t" << safeLog10(eBrem)
											<< "\t" << safeLog10(eDif)
											<< "\t" << safeLog10(eAcc)
											<< "\t" << safeLog10(eAdia) << std::endl;
	

	});


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

	