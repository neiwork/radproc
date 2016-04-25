#include "injection.h"

#include "messages.h"
#include "modelParameters.h"
#include "nonThermalLuminosity.h"
#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\physics.h>




double powerLaw(double E, double Emin, double Emax)
{
	double result = pow(E, (-primaryIndex))*exp(-E / Emax)*exp(-5 * Emin / E);
	return result;
}


double normalization(double z, double t, double Emin, double Emax)
{

	double int_E = RungeKuttaSimple(Emin, Emax, [&Emax, &Emin](double E){
		return E*powerLaw(E, Emin, Emax);
	});  //integra E*Q(E)  entre Emin y Emax

	double vol = pi*P2(jetRadius(rmin, openingAngle))*rmin;  //volumen de la primer celda

	double Q0 = nonThermalLuminosity() / (int_E*vol);  //factor de normalizacion de la inyeccion
	return Q0;
}


void injection(Particle& p, State& st)
{
	show_message(msgStart, Module_electronInjection);

	p.injection.fill([&p, &st](const SpaceIterator& i){
		
		if (i.its[1].canPeek(-1)) 
		{
			return 0.0;
		}
		else //if (z_position = 0) solo inyecto particulas en la posicion 0
		{
			double factor = 0.0, E = i.par.E, z = i.par.R, t = i.par.T;
			//double B = i.par.magneticField;  //VER por que esto no funciona

			double Emin = p.emin();
			double Emax = eEmax(z, magneticField); // 1.6e-12*pow(10.0, p.logEmax);
			double Q0 = normalization(0, 0, Emin, Emax); //VER en principio no depende de z ni t
		
			double total = powerLaw(E, Emin, Emax)*Q0;  //VER calculo el Q0 para no repetir
			return total;
		}

	//	p.injection.set(i, total);
	});

	show_message(msgEnd, Module_electronInjection);
}


