#include "processes.h"




#include "targetFieldM87.h"
#include "write.h"
#include "modelParameters.h"

#include <fluminosities\luminosityHadronic.h>
#include <fluminosities\luminosityPhotoHadronic.h>
#include <fluminosities\luminositySynchrotron.h>
#include <fluminosities\luminosityIC.h>

#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>

#include <fparameters/parameters.h>

#include <fmath\physics.h>
#include <boost/property_tree/ptree.hpp>





double Llab(double Lint, double gamma)
{
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");

	double boost = pow(Dlorentz, 4.0);
	return Lint*boost;
}

/* Takes [emi] =  E^2*[Q(E)] and calculates int(2.0*pi*P2(jetR)*emi dz); 
for [N(E)] = 1/erg, then it just sums over all z and returns erg/s  */




void processes(State& st, const std::string& filename)
{
	
	static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree
	static const double starT = GlobalConfig.get<double>("starT");
	static const double gamma = GlobalConfig.get<double>("Gb");
	static const double Gj = GlobalConfig.get<double>("Gj");
	static const double Lj = GlobalConfig.get<double>("Lj");

	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	double z = GlobalConfig.get<double>("z_peak")*pc;
	

	double Emin = boltzmann*starT / 100.0;


	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);
	
	file << "log(E/eV)"
		<< '\t' << "Synchr"
		<< '\t' << "pg"
		<< '\t' << "IC";

	file << std::endl;
		

	double B = computeMagField(z);
	ParamSpaceValues Bfield(st.photon.ps);
	st.photon.ps.iterate([&](const SpaceIterator &i) {
			
		Bfield.set(i, B);
	});

	double beta_j = beta(Gj);
	double rhoj = Lj / (P2(Gj) * pi*P2(jetRadius(z)) * beta_j*P3(cLight))/protonMass;

	double rhoC = 3.0*solarMass / (4.0*pi*P3(jetRadius(z)))*gamma / protonMass;

	/*ParamSpaceValues rho_j(st.photon.ps);
	st.photon.ps.iterate([&](const SpaceIterator &i) {

		rho_j.set(i, rhoj);
	});*/


	st.photon.ps.iterate([&](const SpaceIterator &i) {


		const double E = i.val(DIM_E);

		double Elab = E*Dlorentz;

		double eSyn = luminositySynchrotron(E, st.proton, i, Bfield);

		double epg = luminosityPhotoHadronic(E, st.proton,
			[&E, &z](double E) {
			return starBlackBody(E, z); }
			, i.coord, Emin, 1.0e4*Emin);
				

		double eIC = luminosityIC(E, st.proton, i,
			[&E, &z](double E) {
			return starBlackBody(E, z); }
		, Emin);

		double epp = luminosityHadronic(E, st.proton, rhoC, i);


		double Lsyn = Llab(eSyn, gamma);
		double Lic = Llab(eIC, gamma);
		double Lpg = Llab(epg, gamma);
		double Lpp = Llab(epp, gamma); 

		double fmtE = log10(Elab / 1.6e-12);
			
			
		file << fmtE
			<< '\t' << safeLog10((Lsyn))
			<< '\t' << safeLog10((Lpg))
			<< '\t' << safeLog10((Lic))
			<< '\t' << safeLog10((Lpp))
			;

		file << std::endl;
				

	});		
	
}
	
	