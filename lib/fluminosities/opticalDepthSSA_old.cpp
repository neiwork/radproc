#include "opticalDepthSSA.h"

#include <fmath\RungeKutta.h>
#include <fparameters\parameters.h>
#include <fmath\interpolation.h>


double fSSA(double x, double E, const Particle& p, const ParamSpaceValues& magf, const SpaceCoord& psc)         //funcion a integrar   x=Ee; L=L(Ega)
{
	//double r = creator.ps.current->par.R;
	//double t = creator.ps.current->par.T;

	double mass = p.mass;
	double Emin = p.emin();
	double Emax = p.emax();

	const double magneticField = magf.get(psc);

	double cte = pow(3.0, 0.5)*P3(electronCharge)*magneticField / (planck*mass*cLight2);

	int n = 100;

	double x_int = pow((Emax / Emin), (1.0 / n));

	////////////// calculo el Psyn
	double Echar = 3.0 * electronCharge*planck*magneticField*P2(x) / (4.0 * pi*P3(mass)*cLight*P2(cLight2));
	double aux = E / Echar;  //aca el aux es el x real

	double Psyn = cte*1.85*pow(aux, (1.0 / 3.0))*exp(-aux);
	//////////////////////////////////
	//  ahora calculo la derivada de N(E)

	double x2 = x*x_int;  // este es el que le sigue
	double dist_x = p.distribution.interpolate({ { 0, x } }, &psc); //creator.dist(x);
	double dist_x2 = p.distribution.interpolate({ { 0, x2 } }, &psc); //creator.dist(x2);  //aca evaluo en x2

	double dx = x*(x_int - 1);

	double deriv = (dist_x2 / P2(x2) - dist_x / P2(x)) / dx;
	double deriv2 = (dist_x2 - dist_x) / (P2(x)*dx) - 2.0*dist_x / P3(x);
	////////////////////////////

	double fSyn = P2(x)*Psyn*deriv;  // funcion a integrar

	return fSyn;
}

double opticalDepthSSA(double E, const Particle& p, const SpaceCoord& psc, const ParamSpaceValues& magf)//(double E, double mass, double Emin, double Emax, const Particle& creator)  //E=Eph
{
	
	double mass = p.mass;
	double Emin = p.emin();
	double Emax = p.emax();

	const double magneticField = magf.get(psc);

	double cte = pow(3.0, 0.5)*P3(electronCharge)*magneticField / (planck*mass*cLight2);
	
	int n = 100;

	double x_int = pow((Emax/Emin),(1.0/n)) ;

	double integral = RungeKuttaSimple(p.emin(), p.emax(), [&](double x) {
		return fSSA(x, E, p, magf, psc);
	});

	double absorptionCoefficient = - P3(planck)*cLight2*integral/(8*pi*P2(E));


	double z = 0.0; //ver!!
	double opacity = absorptionCoefficient*z; // parameters.radius;
	
	//psc[1] ver como obtengo z
	
	return opacity;
}





////double J1(0.0), J2(0.0), J3(0.0);
//double L1 = 0.0;
//
//if (Emin<Emax)
//{
//
//	double x = Emin;
//
//	for (int i = 0; i < n; ++i)
//	{
//		double dx = x*(x_int - 1);
//
//		////////////// calculo el Psyn
//		double Echar = 3.0 * electronCharge*planck*magneticField*P2(x) / (4.0 * pi*P3(mass)*cLight*P2(cLight2));
//		double aux = E / Echar;  //aca el aux es el x real
//
//		double Psyn = cte*1.85*pow(aux, (1.0 / 3.0))*exp(-aux);
//		//////////////////////////////////
//		//  ahora calculo la derivada de N(E)
//
//		double x2 = x*x_int;  // este es el que le sigue
//		double dist_x = p.distribution.interpolate({ { 0, x } }, &psc); //creator.dist(x);
//		double dist_x2 = p.distribution.interpolate({ { 0, x2 } }, &psc); //creator.dist(x2);  //aca evaluo en x2
//
//		double deriv = (dist_x2 / P2(x2) - dist_x / P2(x)) / dx;
//		double deriv2 = (dist_x2 - dist_x) / (P2(x)*dx) - 2.0*dist_x / P3(x);
//		////////////////////////////
//
//		double fSyn = P2(x)*Psyn*deriv;  // funcion a integrar
//
//		L1 = L1 + fSyn*dx;
//
//		//incremento(i,0,n-1,J1,J2,J3,L1);  //if (L1 > 0.0) 
//
//		x = x*x_int;
//
//	}
//}
//
//double integral = L1;// (J1 + 2.0*J2 + 4.0*J3) / 3.0;