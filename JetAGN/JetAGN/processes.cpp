#include "processes.h"

//#include "State.h"

#include "write.h"
#include "messages.h"
#include "luminosityAnisotropicIC.h"
#include <fluminosities\luminositySynchrotron.h>

//#include <fparticle\Particle.h>
//#include <fluminosities/luminosityHadronic.h>
//#include <fluminosities/luminosityPhotoHadronic.h>

#include <omp.h>

#include <fmath\physics.h>
//#include <map>

//class File;

double emiToLumi(const ParamSpace& pps, ParamSpaceValues& psv, int E_position, int t_position);



void processes(State& st)
{
	show_message(msgStart, Module_luminosities);

	ParamSpaceValues Qsyn(st.photon.ps);
	ParamSpaceValues Qic(st.photon.ps);

	std::ofstream file;
	file.open("ntLuminosity_.txt", std::ios::out);
	timestamp_stream(file);

	///////////////////
	//#   pragma omp parallel for \
	//		private(i, eSyn, eIC) \
	//		shared(st, Qsyn, Qic) \
	//		default(none) \
	//		schedule(static, 1) \
	//		num_threads(2)

	st.photon.ps.iterate([&st, &Qsyn, &Qic](const SpaceIterator &i){

		double eSyn = luminositySynchrotron(i.par.E, st.electron); //estos devuelven erg/s/cm^3, integrar!
		double eIC = luminosityAnisotropicIC(i.par.E, st.electron, i.par.R);

		Qsyn.set(i, eSyn);
		Qic.set(i, eIC);
	}, { -1, -1, nR });


	//#pragma omp parallel sections
	//{
	//#pragma omp section
	//	{
	//		Qsyn.fill([&st](const SpaceIterator &i){
	//			return luminositySynchrotron(i.par.E, st.electron); //estos devuelven erg/s/cm^3, integrar!
	//		});
	//	}

	//#pragma omp section
	//	{
	//		Qic.fill([&st](const SpaceIterator &i){
	//			return luminosityAnisotropicIC(i.par.E, st.electron, i.par.R);
	//		});
	//	}
	//}

	const ParamSpace& pps = st.photon.ps;

	timestamp_stream(file);

	file << "log(E/eV)"
		//<< '\t' << "log(t/s)"
		<< '\t' << "log(Lsyn/erg s-1)"
		<< '\t' << "log(Lic/erg s-1)"
		<< std::endl;

//	for (size_t t_ix = 0; t_ix < pps[2].size(); t_ix++) {


		//st.photon.ps.iterate([&st, &Qsyn, &Qic, &file, t_position](const SpaceIterator& i){

		for (size_t E_ix = 0; E_ix < pps[0].size(); E_ix++) {  //estas son las energías del foton, no las del electron

			//E*L(E) = delta^4 E'*L'(E') and E=delta*E'
			//variables primadas ->FF
			//variables sin primar-> lab frame

			//		int E_position = i.its[0].peek; //NO ESTA BIEN ESTO

			int t_ix = nR;
			double E = pps[0][E_ix];
			double t = pps[2][t_ix];  //t es el ultimo valor

			double Elab = E*Dlorentz; //Dlorentz=delta

			double Lsyn = emiToLumi(pps, Qsyn, E_ix, t_ix);
			double Lic = emiToLumi(pps, Qic, E_ix, t_ix);

			double LsynLab = Lsyn*pow(Dlorentz, 4.0);
			double LicLab = Lic*pow(Dlorentz, 4.0);

			double fmtE = log10(Elab / 1.6e-12);

			file << fmtE //<< '\t' << log10(t)
				<< '\t' << safeLog10(LsynLab)
				<< '\t' << safeLog10(LicLab)
				<< std::endl;

		}

	//}
	file.close();

	show_message(msgEnd, Module_luminosities);
}



double emiToLumi(const ParamSpace& pps, ParamSpaceValues& psv, int E_ix, int t_ix) 
{
	double sum = 0.0;
			

	Vector& z = pps.dimensions[1]->values; 

	for (size_t i = 0; i < z.size()-1; ++i) { //no llego al ultimo

//				double dx = x*(x_int - 1);
			double dz = z[i+1] - z[i];

			double jetR = jetRadius(z[i], openingAngle);

			double E = pps[0][E_ix];
			double T = pps[2][t_ix];

			double emissivity = psv.interpolate({ { DIM_E, E }, { DIM_R, z[i] }, { DIM_T, T } });

			double L1 = 2.0*pi*P2(jetR)*emissivity*dz;

			sum = sum + L1;

		}
	
	return sum;
}
