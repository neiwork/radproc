#include "injection.h"

#include "messages.h"
#include <iostream>
#include "modelParameters.h"
#include "nonThermalLuminosity.h"
#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\physics.h>




double powerLaw(double E, double Emin, double Emax)
{
	double result = pow(E, (-parameters.primaryIndex))*exp(-E / Emax)*exp(-5 * Emin / E);
	return result;
}


double normalization(Particle& p, const SpaceCoord& distCoord)
{
	int i_z = distCoord[1];
	double z = p.ps[1][i_z];

	double Emin = p.emin();
	double Emax = eEmax(z, parameters.magneticField);

	double int_E = RungeKuttaSimple(Emin, Emax, [&Emax, &Emin](double E){
		return E*powerLaw(E, Emin, Emax);
	});  //integra E*Q(E)  entre Emin y Emax
	
	double Q0 = dLnt(z) / (int_E);  //factor de normalizacion de la inyeccion
	return Q0;
}


void injection(Particle& p, State& st)
{
	show_message(msgStart, Module_electronInjection);

	//volumen total del jet
	double vol = (pi / 3.0)*(P2(jetRadius(parameters.rmax, parameters.openingAngle))*parameters.rmax
		- P2(jetRadius(parameters.rmin, parameters.openingAngle))*parameters.rmin);

	double z_int = pow((parameters.rmax / parameters.rmin), (1.0 / parameters.nR));


	p.injection.fill([&p, &st, &z_int, &vol](const SpaceIterator& i){
		
		if (i.its[2].canPeek(-1)) 
		{
			return 0.0;
		}
		else //if (t_position = 0) solo inyecto particulas a tiempo 0
		{
			double Emin = p.emin();
			double Emax = eEmax(i.par.R, parameters.magneticField);
			double Q0 = normalization(p, i.coord);
		
			double z = i.par.R;			
			double dz = z*(z_int - 1);
			//volumen de la celda i
			double vol_i = pi*P2(jetRadius(z, parameters.openingAngle))*dz;

			double total = powerLaw(i.par.E, Emin, Emax)*Q0*vol_i / vol;

			return total;
		}

	});

	double Lnt_total = nonThermalLuminosity(parameters.rmin, parameters.rmax);
	
	std::cout << "Lnt total" << '\t' << Lnt_total << std::endl;

	show_message(msgEnd, Module_electronInjection);
}


