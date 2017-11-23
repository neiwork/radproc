#include "radiativeLosses.h"


#include "modelParameters.h"
#include "write.h"


#include <flosses/lossesSyn.h>
#include <flosses/nonThermalLosses.h>
#include <flosses/lossesIC.h>
#include <fparameters/SpaceIterator.h>

#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>


void radiativeLosses(State& st, const std::string& filename)
{
	
//	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
//	static const double starT = GlobalConfig.get<double>("starT");
	//static const double Gamma = GlobalConfig.get<double>("Gamma");
	
	//std::vector<File*> files;
	//OFM out;


	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << "Log(E/eV)" << "\t" << "Log(R/pc)"
		<< "\t" << "Synchr"
		<< "\t" << "IC"
		<< "\t" << "ICAux"
		<< "\t" << "Diff"
		<< "\t" << "Acc"
		//<< "\t" << "Ad" 
		<< std::endl;
	
//	double Emin = 8.17e-7*1.1;
	
	st.electron.ps.iterate([&](const SpaceIterator& i){

		const double magf{ st.magf.get(i) };
		//const double tpf{ st.tpf.get(i) };
		double fmtE = log10(i.val(DIM_E) / 1.6e-12);

		double E = i.val(DIM_E);
		double r = i.val(DIM_R);
//        double theta = i.val(DIM_THETA);
		
		double logR = log10(r/pc);

		double B = magf; // i.par.magneticField;

		double eSyn = lossesSyn(E, B, st.electron) / E;
		
		//double eIC = lossesIC(E, st.electron, 
		//	[&E,&r](double E){
		//	return blackBody(E,r);},   
		//		Emin, 1.0e4*Emin ) / i.val(DIM_E);
		
		//double eIC = lossesIC(E, st.electron, st.tpf, i.coord, Emin, 1.0e4*Emin) /E;
					
				
//		double Reff = 1.0e5;
//		double eDif  = diffusionRate(E, Reff, B);
//		double eAcc = accelerationRate(E, B, accEfficiency);
		//double eAdia = adiabaticLosses(i.val(DIM_E), i.val(DIM_R), vel_lat, Gamma) / i.val(DIM_E);
		
		file << fmtE << "\t" << logR 
							<< "\t" << safeLog10(eSyn) 
							//<< "\t" << safeLog10(eIC)
							//<< "\t" << safeLog10(eDif)
							//<< "\t" << safeLog10(eAcc)
							//<< "\t" << safeLog10(eAdia) 
							<< std::endl;
	

	}, { -1, -1, -1 }); 


	file.close();
	

}