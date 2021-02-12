#include "radiativeLosses.h"


#include "write.h"

#include "ioutil.h"
#include "messages.h"
#include <flosses\lossesSyn.h>
#include <flosses\lossesIC.h>
#include <flosses\lossesHadronics.h>
#include <flosses\lossesPhotoHadronic.h>
#include <flosses\nonThermalLosses.h>

#include <fparameters\SpaceIterator.h>

#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>

void radiativeLosses(State& st, const std::string& folder)
{
	show_message(msgStart, Module_radiativeLosses);

	//static const std::string id = GlobalConfig.get<std::string>("id");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double zInt = GlobalConfig.get<double>("Rdiss");
	static const double density = GlobalConfig.get<double>("density");
	
	std::vector<File*> files;
	OFM out;

	File electronLosses(files, out, getFileName(folder, "electronLosses")); // std::string("electronLosses"));

	//std::ofstream file;
	//file.open(filename.c_str(), std::ios::out);

	//out["electronLosses"]->file << "Log(E/eV)" 
	electronLosses.file << "Log(E/eV)"
		<< "\t" << "Synchr"
		<< "\t" << "IC"
		<< "\t" << "Diff"
		<< "\t" << "Acc"
		<< "\t" << "Ad" << std::endl;

	double vel_lat = cLight*openingAngle;
	
	st.electron.ps.iterate([&](const SpaceIterator& i){

		const double B{ st.magf.get(i) };
		double E = i.val(0); 
		double fmtE = log10(E / 1.6e-12);

		double eSyn = lossesSyn(E, B, st.electron) / E;
		double eIC = lossesIC(E, st.electron,st.tpf, i, st.photon.emin(), st.photon.emax()) / E;  //VER unidades de tpf
			//[&i, &st](double E) {
			//if (E <= st.photon.emax() && E >= st.photon.emin()) {
			//	return st.tpf.interpolate({ { 0, E } });// / (P2(E) *4.0*pi*P2(jetRadius(zInt, openingAngle))*cLight);
			//}
			//else { return 0.0; }; },
			
						
		double eDif  = diffusionRate(E, zInt, B);
		double eAcc = accelerationRate(E, B, accEfficiency);
		double eAdia = adiabaticLosses(E, zInt, vel_lat, Gamma) / E;
		
		//out["electronLosses"]->file << fmtE 
			electronLosses.file << fmtE
							<< "\t" << safeLog10(eSyn) 
							<< "\t" << safeLog10(eIC)
							<< "\t" << safeLog10(eDif)
							<< "\t" << safeLog10(eAcc)
							<< "\t" << safeLog10(eAdia) 
							<< std::endl;
	

	}); //fijo t
	
	//file.close();
	
	File protonLosses(files, out, getFileName(folder, "protonLosses"));

	//out["protonLosses"]->file 
		protonLosses.file << "Log(E/eV)"
						<< "\t" << "Synchr" 
						<< "\t" << "pp"
						<< "\t" << "pgamma"
						<< "\t" << "Diff"
						<< "\t" << "Acc"
						<< "\t" << "Ad" << std::endl;

	st.proton.ps.iterate([&](const SpaceIterator& i){

		double E = i.val(0); 
		double fmtE = log10(E / 1.6e-12);
		double B = st.magf.get(i);

		double pSyn   = lossesSyn(E, B, st.proton) / E;
		double pp     = lossesHadronics(E, density, st.proton)/E;
		double pgamma = lossesPhotoHadronic(E, st.proton, st.tpf, i, st.photon.emin(), st.photon.emax()) / E;
		double pDif   = diffusionRate(E, zInt, B);
		double pAcc   = accelerationRate(E, B, accEfficiency);
		double pAdia = adiabaticLosses(E, zInt, vel_lat, Gamma) / E;

		//out["protonLosses"]->file 
			protonLosses.file << fmtE	<< "\t" << safeLog10(pSyn)
											<< "\t" << safeLog10(pgamma)
											<< "\t" << safeLog10(pp)
											<< "\t" << safeLog10(pDif)
											<< "\t" << safeLog10(pAcc)
											<< "\t" << safeLog10(pAdia) 
											<< std::endl;


	});


}
	