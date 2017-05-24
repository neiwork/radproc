#include "processes.h"

#include "SSC.h"

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




double Llab(double Lint)
{
	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	double boost = pow(Dlorentz, 4.0);
	return Lint*boost;
}

/* Takes [emi] =  E^2*[Q(E)] and calculates int(2.0*pi*P2(jetR)*emi dz); 
for [N(E)] = 1/erg, then it just sums over all z and returns erg/s  */


double emiToLumi(const ParamSpace& pps, ParamSpaceValues& psv, double E, int t_ix)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	double sum = 0.0;

	const double RMIN = pps[DIM_R].first();
	const double RMAX = pps[DIM_R].last();
	const int N_R = pps[DIM_R].size()-1;

	double z_int = pow((RMAX / RMIN), (1.0 / N_R));

	Vector& z = pps.dimensions[1]->values;

	for (size_t i = 0; i < z.size(); ++i) { 

		
		double T = pps[2][t_ix];

		double emissivity = psv.interpolate({ { DIM_E, E }, { DIM_R, z[i] }, { DIM_T, T } });

		//double L1 = 2.0*pi*P2(jetR)*emissivity*dz;
		//la integral de vol la reemplace por una suma sobre todo los z, 
		//ya que N(E) esta multiplicado por el volumen de cada celda

		sum = sum + emissivity;

	}

	return sum;
}


void processes(State& st, const std::string& filename)
{
	static const std::string id = GlobalConfig.get<std::string>("id");

	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	static const double starT = GlobalConfig.get<double>("starT");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	

	show_message(msgStart, Module_luminosities);

	ParamSpaceValues Qsyn(st.photon.ps);
	ParamSpaceValues QicS(st.photon.ps);
	ParamSpaceValues QicAux(st.photon.ps);
	//ParamSpaceValues QicAux2(st.photon.ps);
	
	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	double EphminS = boltzmann*starT / 100.0;
	double EphminAux;

	if (id == "M87") {
		//std::cout << "M87" << std::endl;
		static const double cmbT = GlobalConfig.get<double>("cmbT");
		EphminAux = boltzmann*cmbT / 100.0;
	}
	else {
		static const double IRstarT = GlobalConfig.get<double>("IRstarT");
		EphminAux = boltzmann*IRstarT / 100.0;
	}
		
	

	st.photon.ps.iterate([&](const SpaceIterator &i){

		if (i.coord[0] == 0 && i.coord[1] % 5 == 0){
			std::cout << "height: " << i.coord[1] << std::endl; }
		
		const double E = i.val(DIM_E);
		const double r = i.val(DIM_R);
		double eSyn   = luminositySynchrotron(E, st.electron, i.coord, st.magf); //estos devuelven erg/s, sumar!
		

		//double eIC_Aux2 = luminosityIC(E, st.electron, i.coord, [&E, &r](double E) {
		//	return gxM87(E, r); }, EphminS);

		double eICs, eIC_Aux;

		if (id == "M87") {					
			
			eICs = luminosityIC(E, st.electron, i.coord, [&E, &r](double E) {
				return starBlackBody(E, r); }, EphminS);

			eIC_Aux = luminosityIC(E, st.electron, i.coord, [&E, &r](double E) {
				return cmbBlackBody(E, r); }, EphminAux);
		}
		else {
			eICs = luminosityIC(E, st.electron, i.coord, [&E, &r](double E) {
				return starBlackBody(E, r); }, EphminS);

			eIC_Aux = luminosityIC(E, st.electron, i.coord, [&E, &r](double E) {
				return starIR(E, r); }, EphminAux);
		}

		Qsyn.set(i, eSyn);
		QicS.set(i, eICs);
		QicAux.set(i, eIC_Aux);
		//QicAux2.set(i, eIC_Aux2);

	}, { -1, -1,  (int)st.photon.ps[DIM_R].size() - 1}); //0 });//

	const ParamSpace& pps = st.photon.ps;

	file << "log(E/eV)"
		<< '\t' << "Synchr"
		<< '\t' << "IC"
		<< '\t' << "IC-Aux";
		//<< '\t' << "IC-gx";
	
	file<< std::endl;

	//ParamSpaceValues Qssc(st.photon.ps);
	//SSC(st, Qssc, Qsyn);

	for (int E_ix = 0; E_ix < pps[0].size(); E_ix++) {  

			//E*L(E) = delta^4 E'*L'(E') and E=delta*E'
			//variables primadas ->FF
			//variables sin primar-> lab frame

		//int t_ix = 0;// pps[DIM_R].size() - 1;
		int t_ix = pps[DIM_R].size() - 1;
			double E = pps[0][E_ix]; 
			double t = pps[2][t_ix];  //t es el ultimo valor

			double Elab = E*Dlorentz; //Dlorentz=delta

			double Lsyn   = emiToLumi(pps, Qsyn, E, t_ix);
			double LicS = emiToLumi(pps, QicS, E, t_ix);
			double LicAux = emiToLumi(pps, QicAux, E, t_ix);

			double fmtE = log10(Elab / 1.6e-12);

			file << fmtE //<< '\t' << log10(t)
				<< '\t' << safeLog10(Llab(Lsyn))
				<< '\t' << safeLog10(2.0*Llab(LicS))
				<< '\t' << safeLog10(2.0*Llab(LicAux));
			//el 2.0 en IC es para tener en cuenta la 
			//secci'on eficaz anisotropica en el regimen de thomson
			
			//if (id == "M87") {
			//	double LicAux2 = emiToLumi(pps, QicAux2, E, t_ix);
			//	file << '\t' << safeLog10(2.0*Llab(LicAux2));
			//	}

			file << std::endl;

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



//double dz = z[i]*(z_int - 1);
//volumen de la celda i
//double vol_i = pi*P2(jetRadius(z[i], openingAngle))*dz;;
//double E = pps[0][E_ix];