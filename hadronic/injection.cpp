#include "injection.h"

#include "checkPower.h"
#include "modelParameters.h"

#include <fparameters\parameters.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>

#include <fmath\RungeKutta.h>
#include <fmath\physics.h>

#include <boost/property_tree/ptree.hpp>

#include <iostream>




double dLnt(double z, double Gc, double z_int, double Rs)
{
	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double accEfficiency = GlobalConfig.get<double>("accFraction");
	static const double Gj = GlobalConfig.get<double>("Gj");


	double hj = 1.0;

	double beta_j = beta(Gj);
	double beta_c = beta(Gc);

	double beta_rel = (beta_j - beta_c) / (1.0 - beta_j*beta_c);
	double G_rel = 1.0 / sqrt(1.0 - P2(beta_rel));

	double ramP = Lj / (cLight*pi*P2(jetRadius(z)));

	double cte2 = accEfficiency*ramP*pi*cLight*P2(Rs); // =eta*Lj*So/Sj
																 //double Lsc = cte2*beta_c*Gc*P2(1.0 - beta_c / beta_j);  //hydro
	double Lsc = cte2*beta_rel*G_rel*(G_rel*hj - 1.0) / (P2(Gj)*beta_j);  //hydro

	return Lsc;
}


double powerLaw(double E, double Emin, double Emax)
{
	static const double primaryIndex = GlobalConfig.get<double>("primaryIndex");
	
	double result = pow(E, (-primaryIndex))*exp(-E / Emax)*exp(-5 * Emin / E);
	return result;
}



void injection(Particle& p, State& st)
{

	static const double z = GlobalConfig.get<double>("z_peak")*pc;

	static const double Gb = GlobalConfig.get<double>("Gb");


	double Emin = p.emin();

	double B = computeMagField(z);
	double Rs = jetRadius(z);

	//double E = P2(electronMass*cLight2) / (boltzmann*starT) / Gc[z_ix]; //IC
	//double Q = frad_2(E, z, Gc[z_ix], p);

	double Emax = eEmax(z, B);
	//double Emax = p.emax();

	double int_E = RungeKuttaSimple(Emin, Emax, [&Emax, &Emin](double E) {
		return E*powerLaw(E, Emin, Emax);
	});  //integra E*Q(E)  entre Emin y Emax


	p.injection.fill([&](const SpaceIterator& i) {
		const double E = i.val(DIM_E);

			
		double Q0 = dLnt(z, Gb, z, Rs) / (int_E);  //factor de normalizacion de la inyeccion
		double Q0p = Q0; 


		double total = powerLaw(i.val(DIM_E), Emin, Emax)*Q0p;

		return total;
	});

	double Qinj = computeInjectedPower(st.proton.injection);
}



void distribution(Particle& p, State& st)
{

	static const double Gamma = GlobalConfig.get<double>("Gb");
	double z = GlobalConfig.get<double>("z_peak")*pc;
	

	double tesc = 3.0*z/(Gamma*2.0*cLight); 
	

	p.ps.iterate([&](const SpaceIterator& i) {

		double Ne = p.injection.get(i) * tesc;

		p.distribution.set(i, Ne);


	});
}