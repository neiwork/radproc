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


	double mBH = 1.0e7*solarMass;  //black hole mass
	double rg = mBH*gravitationalConstant / cLight2;

	double z0 = 50.0*rg; //50 * Rschw

	double intRmin = rmin;// z0;
	double intRmax = rmax;// pc;

	//double vol = pi*(P2(jetRadius(intRmax, openingAngle))*intRmax - P2(jetRadius(intRmin, openingAngle))*intRmin) / 3.0;
		
	double vol = pi*P2(jetRadius(rmin, openingAngle))*rmin;  //volumen de la primera celda

	double Q0 = nonThermalLuminosity(intRmin, intRmax) / (int_E*vol);  //factor de normalizacion de la inyeccion
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


