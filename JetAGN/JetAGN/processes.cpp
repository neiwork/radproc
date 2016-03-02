#include "processes.h"

#include "State.h"

#include "write.h"
#include "messages.h"
#include "luminosityAnisotropicIC.h"
#include <fluminosities/luminositySynchrotron.h>
//#include <fluminosities/luminosityHadronic.h>
//#include <fluminosities/luminosityPhotoHadronic.h>
#include <fmath\physics.h>
#include <map>

class File;

double emiToLumi(State& st, ParamSpaceValues& psv, int E_position, int t_position);



void processes(State& st)
{
	show_message(msgStart, Module_luminosities);



	ParamSpaceValues Qsyn(st.photon.ps);
	Qsyn.fill([&st](const SpaceIterator &i){
		double eSyn = luminositySynchrotron(i.par.E, st.electron); //estos devuelven erg/s/cm^3, integrar!
		return eSyn;
	});

	ParamSpaceValues Qic(st.photon.ps);
	Qic.fill([&st](const SpaceIterator &i){
		double eIC = luminosityAnisotropicIC(i.par.E, st.electron, i.par.R);
		return eIC;
	});



	std::ofstream file;
	file.open("ntLuminosity.txt", std::ios::out);

	file << "log(E/eV)"
		<< '\t' << "log(t/s)"
		<< '\t' << "log(Lsyn/erg s-1)"
		<< '\t' << "log(Lic/erg s-1)"
		<< std::endl;

	for (size_t t_position = 0; t_position < st.photon.ps[2].size(); t_position++) {


		//st.photon.ps.iterate([&st, &Qsyn, &Qic, &file, t_position](const SpaceIterator& i){

		for (size_t E_position = 0; E_position < st.photon.ps[0].size(); E_position++) {

			//E*L(E) = delta^4 E'*L'(E') and E=delta*E'
			//variables primadas ->FF
			//variables sin primar-> lab frame

			//		int E_position = i.its[0].peek; //NO ESTA BIEN ESTO

			double E = st.photon.ps[0][E_position];
			double t = st.photon.ps[2][t_position];

			double Elab = E*Dlorentz; //Dlorentz=delta


			double Lsyn = emiToLumi(st, Qsyn, E_position, t_position);
			double Lic = emiToLumi(st, Qic, E_position, t_position);

			double LsynLab = Lsyn*pow(Dlorentz, 4.0);
			double LicLab = Lic*pow(Dlorentz, 4.0);

			double fmtE = log10(Elab / 1.6e-12);

			file << fmtE << '\t' << log10(t)
				<< '\t' << safeLog10(LsynLab)
				<< '\t' << safeLog10(LicLab)
				<< std::endl;

		}

	}
	file.close();
}



	show_message(msgEnd, Module_luminosities);


double emiToLumi(State& st, ParamSpaceValues& psv, int E_position, int t_position) 
{
	double sum;
			

	Vector& z = st.photon.ps.dimensions[1]->values; 

	for (size_t i = 0; i < z.size()-1; ++i) { //no llego al ultimo
	//for (int i = 0; i < n; ++i)
		{
//				double dx = x*(x_int - 1);
			double dz = z[i+1] - z[i];

			double jetR = jetRadius(z[i], openingAngle);

			double E = st.photon.ps[0][E_position];
			double T = st.photon.ps[2][t_position];

			double emissivity = psv.interpolate({ E, z[i], T });

			double L1 = 2.0*pi*P2(jetR)*emissivity*dz;

			sum = sum + L1;

		}
	
	return sum;
}
