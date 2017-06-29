#include "distribution.h"

#include "messages.h"

#include "oneZoneDistribution.h"

#include "injection.h"
#include "losses.h"
//#include "timeDistribution.h"

#include <iostream>
#include <flosses\nonThermalLosses.h>
#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\interpolation.h>
#include <fmath\physics.h>

#include <boost/property_tree/ptree.hpp>

void distribution(Particle& p, State& st, Vector& Gc, Vector& Rc)
{

	show_message(msgStart, Module_electronDistribution);
	double z_0 = p.ps[DIM_R].first();

	/*p.ps.iterate([&](const SpaceIterator& i) {
		const double E = i.val(DIM_E);
		const double z = i.val(DIM_R);
		const double magf = st.magf.get(i);

		oneZoneDistribution(p, st, i, Gc, Rc);
	
		
	}, { 0,-1 });
	*/ 

	static const double Gj = GlobalConfig.get<double>("Gamma");
	p.ps.iterate([&](const SpaceIterator& i) {
		const double E = i.val(DIM_E);
		const double z = i.val(DIM_R);
		const double magf = st.magf.get(i);

		double Emax = (z, magf);
		
		int z_ix = i.coord[DIM_R];
		double bE = losses(E, z, p, st, i, Gc[z_ix]);

		double vc = sqrt(1.0-1.0/P2(Gc[i.coord[DIM_R]]));

		double vjet = sqrt(1.0 - 1.0 / P2(Gj));
		double v_rel = cLight*(vjet - vc);

		double Rs = z / Gc[z_ix]; // Rc[z_ix];

		double tloss = E / bE; //en [s]
		double tesc = 1.0 / escapeRate(Rs, v_rel);  //en [s]

		double dist;

		if (tesc < tloss) {
			dist = p.injection.get(i)*tesc;
		}
		else {
			double integral = RungeKuttaSimple(E, Emax, [&Emax, &p, &z](double Ep) {
			return p.injection.interpolate({ { DIM_E, Ep },{ DIM_R, z } }); });
			dist = integral / bE;
		}

		p.distribution.set(i, dist);

	});


}