#include "processes.h"

#include "targetFields.h"
#include "write.h"
#include "messages.h"
#include "luminosityAnisotropicIC.h"
#include <fluminosities\luminositySynchrotron.h>
#include <fluminosities\luminosityIC.h>

#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>

#include <fparameters/parameters.h>

#include <fmath\physics.h>

#include <boost/property_tree/ptree.hpp>

double tpf(double E, const ParamSpaceValues psv, const SpaceCoord& distCoord)
{
	double result = psv.interpolate({ { 0, E } }, &distCoord);
	return result;
}


double Llab(double Lint)
{
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	return Lint*pow(Dlorentz, 4.0);
}

/* Takes [emi] =  E^2*[Q(E)] and calculates int(2.0*pi*P2(jetR)*emi dz); 
for [N(E)] = 1/erg, then it just sums over all z and returns erg/s  */


double emiToLumi(const ParamSpace& pps, ParamSpaceValues& psv, int E_ix, int t_ix)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle", 0.1);

	double sum = 0.0;

	const double RMIN = pps[DIM_R].first();
	const double RMAX = pps[DIM_R].last();
	const int N_R = pps[DIM_R].size()-1;

	double z_int = pow((RMAX / RMIN), (1.0 / N_R));

	Vector& z = pps.dimensions[1]->values;

	for (size_t i = 0; i < z.size() - 1; ++i) { //no llego al ultimo

		double dz = z[i]*(z_int - 1);
		
		//volumen de la celda i
		double vol_i = pi*P2(jetRadius(z[i], openingAngle))*dz;;

		double E = pps[0][E_ix];
		double T = pps[2][t_ix];

		double emissivity = psv.interpolate({ { DIM_E, E }, { DIM_R, z[i] }, { DIM_T, T } });

		//double L1 = 2.0*pi*P2(jetR)*emissivity*dz;
		//la integral de vol la reemplace por una suma sobre todo los z, 
		//ya que N(E) esta multiplicado por el volumen de cada celda

		sum = sum + emissivity;// *vol_i;

	}

	return sum;
}


void processes(State& st, const std::string& filename)
{
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	static const double starT = GlobalConfig.get<double>("starT");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	show_message(msgStart, Module_luminosities);

	ParamSpaceValues Qsyn(st.photon.ps);
	ParamSpaceValues QicS(st.photon.ps);
	//ParamSpaceValues Qssc(st.photon.ps);

	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);
	//double EphminS(0.0), EphminCMB(0.0);
	//targetPhotonEnergies(EphminS, EphminCMB);

	double EphminS = boltzmann*starT / 100.0;

	st.photon.ps.iterate([&](const SpaceIterator &i){

		const double E = i.val(DIM_E);
		const double r = i.val(DIM_R);
		const double eSyn   = luminositySynchrotron(E, st.electron, i.coord, st.magf); //estos devuelven erg/s, sumar!
		const double eICs = luminosityIC(E, st.electron, i.coord, [&E, &r](double E){
			return starBlackBody(E, r); }, EphminS);

		Qsyn.set(i, eSyn);
		QicS.set(i, eICs);

	}, { -1, -1, (int)st.photon.ps[DIM_R].size()-1 });


	//st.photon.ps.iterate([&](const SpaceIterator &i){

	//	const double E = i.val(DIM_E);
	//	const double r = i.val(DIM_R);

	//	const double eSSC = luminosityIC(E, st.electron, i,
	//		[&Qsyn, &i, &r](double E){return tpf(E, Qsyn, i) / (P2(E) *4.0*pi*P2(jetRadius(r, openingAngle))*cLight); }
	//	, st.photon.emin());

	//		//const double eSSC = luminosityIC(E, st.electron, i.coord, [&](double E){
	//		//	return luminositySynchrotron(E, st.electron, i.coord, st.magf) /
	//		//		(P2(E) *4.0*pi*P2(jetRadius(r, openingAngle))*cLight); }, st.photon.emin());
	//			Qssc.set(i, eSSC);
	//}, { -1, -1, (int)st.photon.ps[DIM_R].size() - 1 });

	const ParamSpace& pps = st.photon.ps;

	file << "log(E/eV)"
		//<< '\t' << "log(t/s)"
		<< '\t' << "log(Lsyn/erg s-1)"
		<< '\t' << "log(Lic/erg s-1)"
		<< '\t' << "log(LiSSC/erg s-1)"
		<< std::endl;

		for (int E_ix = 0; E_ix < pps[0].size(); E_ix++) {  

			//E*L(E) = delta^4 E'*L'(E') and E=delta*E'
			//variables primadas ->FF
			//variables sin primar-> lab frame

			int t_ix = pps[DIM_R].size()-1;
			double E = pps[0][E_ix];
			double t = pps[2][t_ix];  //t es el ultimo valor

			double Elab = E*Dlorentz; //Dlorentz=delta

			double Lsyn   = emiToLumi(pps, Qsyn, E_ix, t_ix);
			double LicS = emiToLumi(pps, QicS, E_ix, t_ix);
			//double Lssc = emiToLumi(pps, Qssc, E_ix, t_ix);

			double fmtE = log10(Elab / 1.6e-12);

			file << fmtE //<< '\t' << log10(t)
				<< '\t' << safeLog10(Llab(Lsyn))
				<< '\t' << safeLog10(Llab(LicS))
				//<< '\t' << safeLog10(Llab(Lssc))
				<< std::endl;

		}

	//}
	file.close();

	show_message(msgEnd, Module_luminosities);
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