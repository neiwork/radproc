#include "oneZone.h"


#include "write.h"
#include "losses.h"
#include "processes.h"
#include "targetFields.h"

#include <flosses\lossesIC.h>
#include <flosses\lossesSyn.h>
#include <fluminosities\luminositySynchrotron.h>
#include <fluminosities\luminosityIC.h>
#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\parameters.h>
#include <fmath\physics.h>
#include <boost/property_tree/ptree.hpp>


double integral(Particle& p, int z_ix, double Eprim)
{
	double suma = 0.0;

	const double Emin = p.ps[DIM_E].first();
	const double Emax = p.ps[DIM_E].last();
	const int nE = p.ps[DIM_E].size() - 1;

	double E_int = pow((Emax / Emin), (1.0 / nE));

	int t_ix = 0;

	p.ps.iterate([&](const SpaceIterator& i) {

		const double E = i.val(DIM_E);
		
		double dE = E*(E_int-1.0);

		if (E >= Eprim)
		{
			suma = suma + p.injection.get(i)*dE;
		}

	}, { -1, z_ix, t_ix });

	return suma;

}




void oneZoneDistribution(Particle& p, State& st)
{
	//static const double vWind = GlobalConfig.get<double>("vWind");
	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double starT = GlobalConfig.get<double>("starT");

	int t_ix = 0;
	//int z_ix = 0;

	double vWind = 2.0e7;
	double Mdot = 1.0e-7*solarMass / yr;


	p.distribution.fill([&](const SpaceIterator& i) {

		const double E = i.val(DIM_E);
		const double r = i.val(DIM_R);
		const double t = i.val(DIM_T);
		const double magf = st.magf.get(i);

		
		double phEmin = boltzmann*starT*1.0e-2;
		double phEmax = boltzmann*starT*1.0e2;

		double lossRate = lossesSyn(E, magf, p)
			+ lossesIC(E, p,
				[&E, &r](double E) {
			return starBlackBody(E, r); }, phEmin, phEmax);


		double Rs = (r*0.1)*sqrt(Mdot*vWind*cLight / (4.0*Lj));

		//double lossRate = losses(E, r, p, st, i); 
		double Tloss = lossRate / E; //[]=s^-1
		double Tesc = Gamma*cLight / (10.0*Rs); //[]=s^-1 t-1' = Gamma* t-1

		if(Tesc > Tloss)
		{
			return  p.injection.get(i)/Tesc;
		}
		else {
			int z_ix = i.coord[1];

			double Q = integral(p, z_ix, E);
			return Q / lossRate;
		}

	}, { -1, -1, t_ix });


}




void processesOneZone(State& st, const std::string& filename)
{
	static const std::string id = GlobalConfig.get<std::string>("id");

	static const double Dlorentz = GlobalConfig.get<double>("Dlorentz");
	static const double starT = GlobalConfig.get<double>("starT");

	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	double EphminS = boltzmann*starT / 100.0;
	
	file << "log(E/eV)"
		<< '\t' << "Synchr"
		<< '\t' << "IC-Aux";

	file << std::endl;


	//for (int E_ix = 0; E_ix < pps[0].size(); E_ix++) {
	st.photon.ps.iterate([&](const SpaceIterator &i) {

		//E*L(E) = delta^4 E'*L'(E') and E=delta*E'
		//variables primadas ->FF
		//variables sin primar-> lab frame

		const double E = i.val(DIM_E);
		const double r = i.val(DIM_R);

		double Elab = E*Dlorentz; //Dlorentz=delta

		double Lsyn = luminositySynchrotron(E, st.electron, i.coord, st.magf);
		double LicAux = luminosityIC(E, st.electron, i.coord, [&E, &r](double E) {
			return starBlackBody(E, r); }, EphminS);

		double LicAux2 = luminosityIC(E, st.electron, i.coord, [&E, &r](double E) {
			return starBlackBody(E, r); }, EphminS);

		double fmtE = log10(Elab / 1.6e-12);

		file << fmtE //<< '\t' << log10(t)
		 	<< '\t' << safeLog10(Llab(Lsyn))
			<< '\t' << safeLog10(2.0*Llab(LicAux))
		    << '\t' << safeLog10(2.0*Llab(LicAux2));
		//el 2.0 en IC es para tener en cuenta la 
		//secci'on eficaz anisotropica en el regimen de thomson

		file << std::endl;

	}, { -1, 0, 0 });

	//}
	file.close();

}

