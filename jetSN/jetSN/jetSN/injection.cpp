#include "injection.h"

#include "messages.h"
#include "nonThermalLuminosity.h"
#include "dynamics.h"
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



void injection(Particle& p, State& st, Vector& Gc, Vector& Rc)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gj = GlobalConfig.get<double>("Gamma");

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
		double beta_c = sqrt(1.0 - 1.0 / P2(Gc[z_ix]));
		double beta_j = sqrt(1.0 - 1.0 / P2(Gj));
		double beta_rel = (beta_j - beta_c) / (1.0 - beta_j*beta_c);
		double G_rel = 1.0 / sqrt(1.0 - P2(beta_rel));

		double B = computeMagField(z, G_rel);
		double Rs = Rc[z_ix];

		double Emax = eEmax(z_0, z, Gc[z_ix], B, Rs);

		double int_E = RungeKuttaSimple(Emin, Emax, [&Emax, &Emin](double E) {
			return E*powerLaw(E, Emin, Emax);
		});  //integra E*Q(E)  entre Emin y Emax


		p.injection.fill([&](const SpaceIterator& i) {
			const double E = i.val(DIM_E);
			const double z = i.val(DIM_R);
			
			double Q0 = dLnt(z, Gc[z_ix], z, Rs) / (int_E);  //factor de normalizacion de la inyeccion
			double Q0p = Q0 / P2(Gc[z_ix]);

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


