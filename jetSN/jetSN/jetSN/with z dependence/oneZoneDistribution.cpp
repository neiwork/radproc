#include "oneZoneDistribution.h"

#include "dynamics.h"
#include "losses.h"

#include <flosses\nonThermalLosses.h>
#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>

#include <fparameters/parameters.h>
#include <fmath/elimiGaussiana.h>

#include <boost/property_tree/ptree.hpp>

void oneZoneDistribution(Particle& p, State& st, const SpaceIterator& si, 
	Vector& Gc, Vector& Rc)
{
	static const double Gj = GlobalConfig.get<double>("Gamma");

	int z_ix = si.coord[DIM_R]; //posicion en la coordenada z
	double z = si.val(DIM_R);

	int nE = p.ps[DIM_E].size();// -1;

	Vector Ee(nE, 0.0);
	Vector Qe(nE, 0.0);
	Vector Ne(nE, 0.0);

	//Vector Ee = particle.energyPoints;
	//Vector injection = particle.injection;
	//int ne = Ee.size() - 1;    //el -1 es para que sea el nEnergyPoints y 

	p.ps.iterate([&](const SpaceIterator& i) {

		//const double z = i.val(DIM_R);

		int e_ix = i.coord[DIM_E]; //posicion en la coordenada z

		Ee[e_ix] = i.val(DIM_E);
		Qe[e_ix] = p.injection.get(i);

	}, { -1, z_ix });

	Matrix a;

	matrixInit(a, Ee.size(), Ee.size() + 1, 0);

	double beta_c = sqrt(1.0 - 1.0 / P2(Gc[z_ix]));
	double beta_j = sqrt(1.0 - 1.0 / P2(Gj));

	double beta_rel = (beta_j - beta_c) / (1.0 - beta_j*beta_c);

	double v_rel = cLight*beta_rel;

	double z_int = p.ps[DIM_R].first();
	double Rs = z / Gc[z_ix]; // Rc[z_ix];
	//double Rs = Rc[z_ix] / Gc[z_ix]; 
	
	nE = nE - 1;

	for (int i = 0; i <= nE; ++i) {     //ver como acomodo los limites

		if (i == nE) {
			a[i][i] = 1;    //esta es la condicion inicial: N(Emax)=0
			a[i][nE + 1] = 0;
		}
		else {
			double h = Ee[i + 1] - Ee[i];

			for (int j = 0; j <= nE + 1; ++j) {
				double tesc = 1.0 / escapeRate(Rs, v_rel);  //en [s]  
				double frac = 0.5*h / 1.0e20;
				if (j == i) {                               
					double bE = losses(Ee[i], z, p, st, si, Gc[z_ix]);
					
					//a[i][j] = losses(Ee[i], particle, state) + 0.5*h / escapeTime(Ee[i], particle);
					a[i][j] = bE + 0.5*h / tesc;
				}
				else if (j == i + 1) {
					double bE = losses(Ee[i + 1], z, p, st, si, Gc[z_ix]);
					//a[i][j] = -losses(Ee[i + 1], particle, state) + 0.5*h / escapeTime(Ee[i], particle);
					a[i][j] = -bE + 0.5*h / tesc;
				}
				else if (j == nE + 1) {

					a[i][j] = 0.5*h*(Qe[i] + Qe[i + 1]);
				}
				else {
					a[i][j] = 0;
				}
			}
		}
	}


	elimiGaussiana(nE+1, a, Ne); //ne+1

	p.ps.iterate([&](const SpaceIterator& i) {

		int e_ix = i.coord[DIM_E];

		p.distribution.set(i, Ne[e_ix]);

	}, { -1, z_ix });

}