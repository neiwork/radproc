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

double Llab(double Lint)
{
	static const double Dlorentz = GCFG.get<double>("Dlorentz");
	return Lint*pow(Dlorentz, 4.0);
}

/* Takes [emi] =  E^2*[Q(E)] and calculates int(2.0*pi*P2(jetR)*emi dz); 
for [N(E)] = 1/erg, then it just sums over all z and returns erg/s  */


double emiToLumi(const ParamSpace& pps, ParamSpaceValues& psv, int E_ix, int t_ix)
{
	static const double openingAngle = GCFG.get<double>("openingAngle", 0.1);

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

		sum = sum + emissivity*vol_i;

	}

	return sum;
}


void processes(State& st, const std::string& filename)
{
	static const double Dlorentz = GCFG.get<double>("Dlorentz");

	show_message(msgStart, Module_luminosities);

	ParamSpaceValues Qsyn(st.photon.ps);
	ParamSpaceValues Qic(st.photon.ps);
	ParamSpaceValues QicCMB(st.photon.ps);

	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);
	timestamp_stream(file);

	double EphminS(0.0), EphminCMB(0.0);
	targetPhotonEnergies(EphminS, EphminCMB);

	st.photon.ps.iterate([&](const SpaceIterator &i){

		const double E = i.val(DIM_E);
		const double eSyn   = luminositySynchrotron(E, st.electron, i.coord, st.magf); //estos devuelven erg/s/cm^3, integrar!
		const double eIC    = luminosityAnisotropicIC(E, st.electron, i.coord, EphminS);
		const double eICcmb = luminosityIC( E, st.electron, i.coord, cmbBlackBody, EphminCMB);

		Qsyn.set(i, eSyn);
		Qic.set(i, eIC);
		QicCMB.set(i, eICcmb);

	}, { -1, -1, (int)st.photon.ps[DIM_R].size()-1 });



	const ParamSpace& pps = st.photon.ps;

	timestamp_stream(file);

	file << "log(E/eV)"
		//<< '\t' << "log(t/s)"
		<< '\t' << "log(Lsyn/erg s-1)"
		<< '\t' << "log(Lic/erg s-1)"
		<< '\t' << "log(LicCMB/erg s-1)"
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
			double Lic    = emiToLumi(pps, Qic, E_ix, t_ix);
			double LicCMB = emiToLumi(pps, QicCMB, E_ix, t_ix);


			double fmtE = log10(Elab / 1.6e-12);

			file << fmtE //<< '\t' << log10(t)
				<< '\t' << safeLog10(Llab(Lsyn))
				<< '\t' << safeLog10(Llab(Lic))
				<< '\t' << safeLog10(Llab(LicCMB))
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