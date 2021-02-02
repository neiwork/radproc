#include "opticalDepthSSA.h"

#include <fmath/RungeKutta.h>
#include <fparameters/parameters.h>
#include <fmath/interpolation.h>


double fSSA(double x, double E, const Particle& p, const double magneticField, const SpaceCoord& psc)         //funcion a integrar   x=Ee; L=L(Ega)
{
	double mass = p.mass;	
	double cte = pow(3.0, 0.5)*P3(electronCharge)*magneticField / (planck*mass*cLight2);

	////////////// calculo el Psyn
	double Echar = 3.0 * electronCharge*planck*magneticField*P2(x) / (4.0 * pi*P3(mass)*cLight*P2(cLight2));
	double aux = E / Echar;  //aca el aux es el x real

	double Psyn = cte*1.85*pow(aux, (1.0 / 3.0))*exp(-aux);
	//////////////////////////////////
	
//****  calculo la derivada de N(E)
	
	double Emin = p.emin();
	double Emax = p.emax();
	
	int n = 30; //opt.samples_x; //misma cantidad de puntos que el RK

	double x_int = pow((Emax / Emin), (1.0 / n));

	double x1 = x/sqrt(x_int);
	double x2 = x*sqrt(x_int);  // este es el que le sigue
	//double dist_x = p.distribution.interpolate({ { 0, x } }, &psc);
	//double dist_x2 = p.distribution.interpolate({ { 0, x2 } }, &psc); 
	
	double dist_x, dist_x1, dist_x2;
	if (x < Emin || x> Emax){
		dist_x = 0.0;
	}
	else{
		dist_x = p.distribution.interpolate({ { 0, x } }, &psc);
	}
	
	if (x1 < Emin || x1 > Emax) {
		dist_x1 = 0.0;
	}
	else {
		dist_x1 = p.distribution.interpolate({{0,x1}},&psc);
	}

	if (x2 < Emin || x2> Emax) {
		dist_x2 = 0.0;
	}
	else{
		dist_x2 = p.distribution.interpolate({ { 0, x2 } }, &psc);
	}

	double dx = x1*(x_int - 1);

	double deriv = (dist_x2 / P2(x2) - dist_x1 / P2(x1)) / dx;
	//double deriv2 = (dist_x2 - dist_x) / (P2(x)*dx) - 2.0*dist_x / P3(x);
	////////////////////////////

	return P2(x)*Psyn*deriv;  // funcion a integrar
}

double fSSA2(double x, double E, const Particle& p, const double magf, 
				const SpaceCoord& psc)
{
	double cte = pow(3.0, 0.5)*P3(electronCharge)*magf/(planck*p.mass*cLight2);

	////////////// calculo el Psyn
	double Echar = 3.0*electronCharge*planck*magf*P2(x)/(4.0*pi*P3(p.mass)*cLight*P2(cLight2));
	double aux = E/Echar;  //aca el aux es el x real

	double Psyn = (magf > 0.0) ? cte*1.85*pow(aux, (1.0 / 3.0))*exp(-aux) : 0.0;
	//////////////////////////////////
	
//****  calculo la derivada de N(E)
	
	double Emin = p.emin();
	double Emax = p.emax();
	int n = 30; //opt.samples_x; //misma cantidad de puntos que el RK
	double x_int = pow((p,Emax/Emin),1.0/n);

	double x1 = x/sqrt(x_int);
	double x2 = x*sqrt(x_int);
	
	double dist_x = (x > Emin && x < Emax) ? p.distribution.interpolate({{0,x}},&psc) : 0.0;
	double dist_x1 = (x1 > Emin && x1 < Emax) ? p.distribution.interpolate({{0,x1}},&psc) : 0.0;
	double dist_x2 = (x2 > Emin && x2 < Emax) ? p.distribution.interpolate({{0,x2}},&psc) : 0.0;

	double dx = x1*(x_int - 1.0);

	double deriv = (dist_x2/P2(x2) - dist_x1/P2(x1)) / dx;
	return x*x*Psyn*deriv;  // funcion a integrar
}

//esta rutina esta implementada en el proyecto porque depende de las coordenadas
/*double opticalDepthSSA(double E, const Particle& p, const SpaceCoord& psc, const double magneticField)//(double E, double mass, double Emin, double Emax, const Particle& creator)  //E=Eph
{
	

	double mass = p.mass;

	//const double magneticField = magf.get(psc);

	double cte = pow(3.0, 0.5)*P3(electronCharge)*magneticField / (planck*mass*cLight2);

	double integral = RungeKuttaSimple(p.emin(), p.emax(), [&](double x) {
		return fSSA(x, E, p, magneticField, psc);
	});

	double absorptionCoefficient = - P3(planck)*cLight2*integral/(8*pi*P2(E));


	double z = 0.0; //ver!!
	double opacity = absorptionCoefficient*z; // parameters.radius;
	
	return opacity;
}*/

