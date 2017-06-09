#include "injection.h"

#include "messages.h"
#include "modelParameters.h"
#include "nonThermalLuminosity.h"

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



void injection(Particle& p, State& st, Vector& Gc)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	double z_0 = p.ps[DIM_R].first();

	show_message(msgStart, Module_electronInjection);

	const double RMIN = p.ps[DIM_R].first();
	const double RMAX = p.ps[DIM_R].last();
	const int N_R = p.ps[DIM_R].size()-1;

	//double Lnt_total = nonThermalLuminosity(p, Gc);

	//volumen total del jet
	//double vol = (pi / 3.0)*(P2(jetRadius(RMAX, openingAngle))*RMAX
	//	- P2(jetRadius(RMIN, openingAngle))*RMIN);

	double R_int = pow((RMAX / RMIN), (1.0 / N_R));
	
	double Emin = p.emin();

	double z = RMIN;

	for (int z_ix = 0; z_ix < N_R; z_ix++) {

		double dz = z*(R_int - 1);

		//normalizacion afuera del iterate sobre E
		double B = computeMagField(z);
		double Emax = eEmax(z, B);

		double int_E = RungeKuttaSimple(Emin, Emax, [&Emax, &Emin](double E) {
			return E*powerLaw(E, Emin, Emax);
		});  //integra E*Q(E)  entre Emin y Emax

		double Q0 = dLnt(z, Gc[z_ix], z_0) / (int_E);  //factor de normalizacion de la inyeccion
		double Q0p = Q0 / P2(Gc[z_ix]); 
		//////

		p.injection.fill([&](const SpaceIterator& i) {
			const double E = i.val(DIM_E);
			const double z = i.val(DIM_R);

			//double z = i.val(DIM_R);			
			//double dz = z*(z_int - 1);
			//volumen de la celda i
			double vol_i = pi*P2(jetRadius(z, openingAngle))*dz;

			double total = powerLaw(i.val(DIM_E), Emin, Emax)*Q0p;// *vol_i;// / vol;

			return total;
		}, { -1, z_ix });
	
	}
	show_message(msgEnd, Module_electronInjection);
}


