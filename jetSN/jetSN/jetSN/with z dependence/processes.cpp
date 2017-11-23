#include "processes.h"


#include "checkPower.h"
//#include "SSC.h"

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





double Llab(double Lint, double gamma)
{
	double Dlorentz = computeDlorentz(gamma);
	double boost = pow(Dlorentz, 4.0);
	return Lint*boost;
}

/* Takes [emi] =  E^2*[Q(E)] and calculates int(2.0*pi*P2(jetR)*emi dz); 
for [N(E)] = 1/erg, then it just sums over all z and returns erg/s  */


void synL(State& st, ParamSpaceValues& Qsyn)
{


	st.photon.ps.iterate([&](const SpaceIterator &i) {

		const double E = i.val(DIM_E);
		double eSyn = luminositySynchrotron(E, st.electron, i.coord, st.magf); //estos devuelven erg/s, sumar!

		Qsyn.set(i, eSyn);
	});
}

void processes(State& st, const std::string& filename, Vector& Gc, Vector& tobs)
{
	
	static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree
	static const double starT = GlobalConfig.get<double>("starT");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	
	ParamSpaceValues Qsyn(st.photon.ps);
	synL(st, Qsyn);

	show_message(msgStart, Module_luminosities);

	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	std::ofstream file2;
	file2.open(filename + "_peak_evol.txt", std::ios::out);

	double EphminS = boltzmann*starT / 100.0;

	static const double IRstarT = GlobalConfig.get<double>("IRstarT");
	double EphminAux = boltzmann*IRstarT / 100.0;

	const double z0 = st.electron.ps[DIM_R].first();
	double t0 = z0 / cLight;
	

	file << "log(E/eV)"
		<< '\t' << "z/pc"
		<< '\t' << "tobs/yr"
		<< '\t' << "Synchr"
		<< '\t' << "IC"
		<< '\t' << "IC - IR"
	    << '\t' << "SSC";

	file << std::endl;

	const int N_R = st.photon.ps[DIM_R].size() - 1;
	
	for (int z_ix = 0; z_ix < N_R; z_ix++) {

		double t = tobs[z_ix];

		double Ltot = 0.0;

		double gamma = Gc[z_ix];
		double Dlorentz = computeDlorentz(gamma);
		double r = st.photon.distribution.ps[DIM_R][z_ix];

		double E_kn = P2(electronMass*cLight2) / (boltzmann*IRstarT) / gamma;

		SpaceCoord j(st.photon.ps);
		j[DIM_R] = z_ix;

		double L_kn = luminositySynchrotron(E_kn, st.electron, j, st.magf)
			 + luminosityIC(E_kn, st.electron, j, [&E_kn, &r, &gamma](double E_kn) {
			return starIR(E_kn, r, gamma); }, EphminAux);

		file2 << t / yr << '\t' << safeLog10(Llab(L_kn, gamma)) << std::endl;

		st.photon.ps.iterate([&](const SpaceIterator &i) {

		//	if (i.coord[DIM_E] == 0 && i.coord[DIM_R] % 5 == 0) {
		//		std::cout << "z: " << i.coord[DIM_R] << std::endl;
		//	}
			

			const double E = i.val(DIM_E);
			//const double r = i.val(DIM_R);

			double Elab = E*Dlorentz; 

			double eSyn = Qsyn.get(i);// luminositySynchrotron(E, st.electron, i.coord, st.magf); //estos devuelven erg/s, sumar!

			double eICs = luminosityIC(E, st.electron, i.coord, [&E, &r, &gamma](double E) {
				return starBlackBody(E, r, gamma); }, EphminS);

			double eIC_aux = luminosityIC(E, st.electron, i.coord, [&E, &r, &gamma](double E) {
				return starIR(E, r, gamma); }, EphminAux);

			double eSSC = luminosityIC(E, st.electron, i,
				[&Qsyn, &i, &r](double E) {return Qsyn.interpolate({ {DIM_E, E },{ DIM_R, r } }) / (P2(E) *4.0*pi*P2(jetRadius(r, openingAngle))*cLight); }
			, st.photon.emin());

			double Ltot1 = eSyn + eICs + eIC_aux + eSSC;
			Ltot = Ltot + Ltot1;

			double Lsyn = Llab(eSyn, gamma);
			double LicS = Llab(eICs, gamma);
			double Lic_aux = Llab(eIC_aux, gamma);
			double Lssc = Llab(eSSC, gamma);

			double fmtE = log10(Elab / 1.6e-12);
			
			//double tlab = (r / cLight - t0) / yr;
			
			
			file << fmtE
				<< '\t' << r / pc
				<< '\t' << t/yr
				<< '\t' << safeLog10((Lsyn))
				<< '\t' << safeLog10((LicS))
				<< '\t' << safeLog10((Lic_aux))
				<< '\t' << safeLog10((Lssc))
				<< '\t' << safeLog10(Llab(Ltot1, gamma))
				;

			file << std::endl;

			

		}, { -1, z_ix });

		double Qinj = computeInjectedPower(st.electron.injection, z_ix);

		std::cout <<  st.photon.ps[DIM_R][z_ix] / pc << '\t' <<
			//"injected power:" << '\t' << Qinj << '\t' 
			//<< "One zone power" << '\t' <<  Ltot
			"oneZone/InjPower:" << '\t' << Ltot/Qinj 
			<< std::endl;

	}

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