#include "injection.h"

#include "messages.h"
#include "modelParameters.h"
#include "nonThermalLuminosity.h"

#include <fparameters\parameters.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>

#include <fmath\RungeKutta.h>
#include <fmath\physics.h>

#include <iostream>


double powerLaw(double E, double Emin, double Emax)
{
	static const double primaryIndex = GCFG.get<double>("primaryIndex", 2.0);
	
	double result = pow(E, (-primaryIndex))*exp(-E / Emax)*exp(-5 * Emin / E);
	return result;
}


double normalization(Particle& p, double z, double magf)
{
	//int i_z = distCoord[1];
	//double z = p.ps[1][i_z];

	double Emin = p.emin();
	double Emax = eEmax(z, magf);

	double int_E = RungeKuttaSimple(Emin, Emax, [&Emax, &Emin](double E){
		return E*powerLaw(E, Emin, Emax);
	});  //integra E*Q(E)  entre Emin y Emax
	
	double Q0 = dLnt(z) / (int_E);  //factor de normalizacion de la inyeccion
	return Q0;
}


void injection(Particle& p, State& st)
{
	static const double openingAngle = GCFG.get<double>("openingAngle", 0.1);

	show_message(msgStart, Module_electronInjection);

	const double RMIN = p.ps[DIM_R].first();
	const double RMAX = p.ps[DIM_R].last();
	const int N_R = p.ps[DIM_R].size()-1;

	//volumen total del jet
	double vol = (pi / 3.0)*(P2(jetRadius(RMAX, openingAngle))*RMAX
		- P2(jetRadius(RMIN, openingAngle))*RMIN);

	double z_int = pow((RMAX / RMIN), (1.0 / N_R));


	p.injection.fill([&p, &st, &z_int, &vol](const SpaceIterator& i){
		const double magf{ st.magf.get(i) };
		const double r{ i.val(DIM_R) };
		if (i.its[2].canPeek(-1)) 
		{
			return 0.0;
		}
		else //if (t_position = 0) solo inyecto particulas a tiempo 0
		{
			double Emin = p.emin();
			double Emax = eEmax(r,magf);
			double Q0 = normalization(p,r,magf);
		
			double z = i.val(DIM_R);			
			double dz = z*(z_int - 1);
			//volumen de la celda i
			double vol_i = pi*P2(jetRadius(z, openingAngle))*dz;

			double total = powerLaw(i.val(DIM_E), Emin, Emax)*Q0*vol_i / vol;

			return total;
		}

	});

	double Lnt_total = nonThermalLuminosity(RMIN, RMAX);
	
	std::cout << "Lnt total" << '\t' << Lnt_total << std::endl;

	show_message(msgEnd, Module_electronInjection);
}


