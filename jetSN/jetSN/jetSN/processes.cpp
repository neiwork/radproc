#include "processes.h"

#include "SSC.h"

#include "targetFields.h"
#include "write.h"
#include "messages.h"
//#include "luminosityAnisotropicIC.h"
#include <fluminosities\luminositySynchrotron.h>
#include <fluminosities\luminosityIC.h>

#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>

#include <fparameters/parameters.h>

#include <fmath\physics.h>

#include <boost/property_tree/ptree.hpp>




double Llab(double Lint)
{
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	double boost = pow(Dlorentz, 4.0);
	return Lint*boost;
}

/* Takes [emi] =  E^2*[Q(E)] and calculates int(2.0*pi*P2(jetR)*emi dz); 
for [N(E)] = 1/erg, then it just sums over all z and returns erg/s  */






void synL(State& st, ParamSpaceValues& Qsyn)
{
	static const double z_int = GlobalConfig.get<double>("z_int")*pc;


	st.photon.ps.iterate([&](const SpaceIterator &i) {

		const double E = i.val(DIM_E);
		double eSyn = luminositySynchrotron(E, st.electron, i.coord, st.magf); //estos devuelven erg/s, sumar!

		Qsyn.set(i, eSyn);
	});
}

void processes(State& st, const std::string& filename)
{
	
	//static const double Bfield = GlobalConfig.get<double>("Bfield");
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	static const double starT = GlobalConfig.get<double>("starT");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double z_int = GlobalConfig.get<double>("z_int")*pc;


	ParamSpaceValues Qsyn(st.photon.ps);
	synL(st, Qsyn);

	show_message(msgStart, Module_luminosities);

	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	double EphminS = boltzmann*starT / 100.0;

	static const double IRstarT = GlobalConfig.get<double>("IRstarT");
	double EphminAux = boltzmann*IRstarT / 100.0;


	file << "log(E/eV)"
		<< '\t' << "Synchr"
		<< '\t' << "IC"
		<< '\t' << "IC-Aux"
	    << '\t' << "SSC";

	file << std::endl;

	st.photon.ps.iterate([&](const SpaceIterator &i) {

		if (i.coord[0] % 5 == 0){
			std::cout << "energy: " << i.coord[0] << std::endl; }

		const double E = i.val(DIM_E);
		double Elab = E*Dlorentz; //Dlorentz=delta

		const double r = z_int;
		double eSyn = Qsyn.get(i);// luminositySynchrotron(E, st.electron, i.coord, st.magf); //estos devuelven erg/s, sumar!

		double eICs = luminosityIC(E, st.electron, i.coord, [&E, &r](double E) {
			return starBlackBody(E, r); }, EphminS);

		double eIC_aux = luminosityIC(E, st.electron, i.coord, [&E, &r](double E) {
			return starIR(E, r); }, EphminAux);

		double eSSC = luminosityIC(E, st.electron, i, 
		[&Qsyn, &i, &r](double E) {return Qsyn.interpolate({ {DIM_E, E } }) / (P2(E) *4.0*pi*P2(jetRadius(r, openingAngle))*cLight); }
		, st.photon.emin());


		double Lsyn = Llab(eSyn);
		double LicS = Llab(eICs);
		double Lic_aux = Llab(eIC_aux);
		double Lssc = Llab(eSSC);

		double fmtE = log10(Elab / 1.6e-12);

		file << fmtE //<< '\t' << log10(t)
			<< '\t' << safeLog10((Lsyn))
			<< '\t' << safeLog10((LicS))
			<< '\t' << safeLog10((Lic_aux))
			<< '\t' << safeLog10((Lssc))
			;


		file << std::endl;

	});

}
	
	





///////////////////
//#   pragma omp parallel for \
	//		private(i, eSyn, eIC) \
	//		shared(st, Qsyn, Qic) \
	//		default(none) \
	//		schedule(static, 1) \
	//		num_threads(2)

//#pragma omp parallel sections
//{
//#pragma omp section
//	{
//		Qsyn.fill([&st](const SpaceIterator &i){
//			return luminositySynchrotron(i.val(DIM_E), st.electron); //estos devuelven erg/s/cm^3, integrar!
//		});
//	}

//#pragma omp section
//	{
//		Qic.fill([&st](const SpaceIterator &i){
//			return luminosityAnisotropicIC(i.val(DIM_E), st.electron, i.val(DIM_R));
//		});
//	}
//}



//double dz = z[i]*(z_int - 1);
//volumen de la celda i
//double vol_i = pi*P2(jetRadius(z[i], openingAngle))*dz;;
//double E = pps[0][E_ix];