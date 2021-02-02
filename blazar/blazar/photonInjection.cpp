#include "photonInjection.h"

#include "modelParameters.h"
#include "ebl_absorption.h"

#include <fparameters\parameters.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>

#include <boost/property_tree/ptree.hpp>

#include <finjection\neutronPgamma.h>
#include <finjection\pgammaPionInj.h>
#include <fluminosities\luminositySynchrotron.h>
#include <fluminosities\luminosityIC.h>
#include <fluminosities\luminosityPhotoHadronic.h>
#include <fluminosities\luminosityHadronic.h>
#include <fparameters/parameters.h>

#include <string>
#include <iostream>
#include <fstream>


//double tpf(double E, const ParamSpaceValues psv, const SpaceCoord& distCoord)
//double tpf(double E, Particle& photon, const SpaceCoord& distCoord)
//{
//	const ParamSpaceValues psv = photon.injection;
//	double result = 0.0;
//	double emin = photon.emin();
//
//	if (E >= emin && E <= photon.emax()) {
//		result = psv.interpolate({ { 0, E } }, &distCoord);
//	}
//	return result;
//}


double Llab(double Lint, double Dlorentz)
{
	//double Dlorentz = computeDlorentz(gamma);
	double boost = pow(Dlorentz, 4.0);
	return Lint*boost;
}



//void photonDistribution(Particle& p, State& st)
//{
void processes(State& st, const std::string& filename)
{

	//static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree
	//static const double Gamma = GlobalConfig.get<double>("Gamma");
	//static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	static const double density = GlobalConfig.get<double>("density");

	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	double magf = st.magf.get({ 0 });
	//double Emin = p.emin();
	//double Emax = p.emax();//es la del foton eEmax(zInt, B);


	Vector E_TeV(20, 0.0);
	Vector tau02(20, 0.0);
	Vector tau1(20, 0.0);

	absorption(E_TeV, tau02, tau1);
	
	file << "log(E/eV)"
		<< '\t' << "Synchr"
		<< '\t' << "SSC";

	file << std::endl;
		
	//double Dlorentz = computeDlorentz(gamma);


	//double B = st.magf.get({ 0 });
	//double Emin = p.emin();
	//double Emax = eEmax(Rdiss, B);

	//volumen 
	//double vol = pi*P2(jetRadius(Rdiss, openingAngle))*Rdiss;
	//no es necesario el volumen, [Q]= 1/ erg s

	//const ParamSpaceValues psv = st.photon.injection;
	
	//ParamSpaceValues Qsyn(st.photon.ps);
	//synL(st, Qsyn);

	st.photon.distribution.ps.iterate([&](const SpaceIterator& i){

		double E = i.val(DIM_E);
		double Elab = E*Dlorentz;
		
		pgammaPionInj(E, st.proton, st.tpf, i, st.photon.emin(), st.photon.emax());
		neutronPgamma(E, st.proton, st.proton, st.tpf, i, st.photon.emin(), st.photon.emax());

		double eSyn2 = luminositySynchrotron2(E, st.electron, i, st.magf.get(i));
		//double eSyn = luminositySynchrotron(E, st.electron, i, st.magf);

		double eSSC = luminosityIC_old(E, st.electron, i, 
			[&i, &st](double E) {
			if (E <= st.photon.emax() && E >= st.photon.emin()) {
				return st.tpf.interpolate({ { DIM_E, E } });// / (P2(E) *4.0*pi*P2(jetRadius(zInt, openingAngle))*cLight);
			}
			else { return 0.0; }; }
			, st.photon.emin(), st.photon.emax());

		double pPP = luminosityHadronic(E, st.proton, density, i);

		double pPG = luminosityPhotoHadronic(E, st.proton, st.tpf, i, st.photon.emin(), st.photon.emax());

		double Ltot_int = eSyn2 + eSSC;


		double Lsyn = Llab(eSyn2, Dlorentz);
		double Lssc = Llab(eSSC, Dlorentz);
		double Ltot = Llab(Ltot_int, Dlorentz);

		double fmtE = log10(Elab / 1.6e-12);

		st.photon.injection.set(i,Lsyn+Lssc);

		double e_tau02 = 1.0;
		double e_tau1 = 1.0;

		double ElabTeV = Elab / 1.6;
		if (ElabTeV >= 0.1) {
			e_tau02 = interpol(ElabTeV, E_TeV, tau02, E_TeV.size(), 0);
			e_tau1 = interpol(ElabTeV, E_TeV, tau1, E_TeV.size(), 0);
		}

		file << fmtE
			//	<< '\t' << r / pc
			//	<< '\t' << t/yr
			<< '\t' << safeLog10((Lsyn))
			<< '\t' << safeLog10((Lssc))
			<< '\t' << safeLog10(Ltot)
			<< '\t' << safeLog10(Ltot*e_tau02)
			<< '\t' << safeLog10(Ltot*e_tau1)
			;

		file << std::endl;

	});


}