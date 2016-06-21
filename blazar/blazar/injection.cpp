#include "injection.h"


#include "modelParameters.h"

#include <fparameters\parameters.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>

#include <fmath\RungeKutta.h>
#include <fmath\physics.h>

#include <boost/property_tree/ptree.hpp>

#include <iostream>



double powerLaw(double E, double Emin, double Emax)
{
	static const double primaryIndex = GlobalConfig.get<double>("primaryIndex");
	
	double result = pow(E, (-primaryIndex))*exp(-E / Emax)*exp(-5 * Emin / E);
	return result;
}


double normalization(Particle& p, double z, double magf)
{
	static const double Gamma = GlobalConfig.get<double>("Gamma");

	double eta = GlobalConfig.get<double>("accEfficiency");
	double Lj = GlobalConfig.get<double>("Lj");

	double Emin = p.emin();
	double Emax = eEmax(z, magf);

	double int_E = RungeKuttaSimple(Emin, Emax, [&Emax, &Emin](double E){
		return E*powerLaw(E, Emin, Emax);
	});  //integra E*Q(E)  entre Emin y Emax
	
	double Q0 = eta*Lj / (int_E);  //factor de normalizacion de la inyeccion
	return Q0 / Gamma;
	//N'(E')dE' = N(E)dE  ==> E'_nt = int( E'N(E')dE') = E_nt/Gamma;
}


void injection(Particle& p, State& st)
{
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	static const double Rdiss = GlobalConfig.get<double>("Rdiss");

	double B = st.magf.get({ 0 });
	double Emin = p.emin();
	double Emax = eEmax(Rdiss, B);

	double Q0 = normalization(p, Rdiss, B);
	
	//volumen 
	double vol = pi*P2(jetRadius(Rdiss, openingAngle))*Rdiss;


	p.injection.fill([&](const SpaceIterator& i){
		
		double total = powerLaw(i.val(DIM_E), Emin, Emax)*Q0;// *vol;

		return total;
	});


}


