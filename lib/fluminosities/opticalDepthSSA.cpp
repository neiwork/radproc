#include "opticalDepthSSA.h"

#include <fparameters\parameters.h>
//#include <fmath\interpolation.h>




double opticalDepthSSA(double E, double mass, double Emin, double Emax, const Particle& creator)  //E=Eph
{
	
	double cte	= pow(3.0,0.5)*P3(electronCharge)*magneticField/(planck*mass*cLight2);
	
	int n = 100;

	double x_int = pow((Emax/Emin),(1.0/n)) ;

	//double J1(0.0), J2(0.0), J3(0.0);
	double L1 = 0.0;

	if(Emin<Emax)
	{
		double x = Emin;
		
		for (int i=0; i < n; ++i)
		{
			double dx = x*(x_int-1);

////////////// calculo el Psyn
			double Echar = 3*electronCharge*planck*magneticField*P2(x)/(4*pi*P3(mass)*cLight*P2(cLight2));
			double aux = E/Echar;  //aca el aux es el x real

			double Psyn = cte*1.85*pow(aux,(1.0/3.0))*exp(-aux);  
//////////////////////////////////
//  ahora calculo la derivada de N(E)

			double x2 = x*x_int;  // este es el que le sigue
			double dist_x  = creator.dist(x);
			double dist_x2 = creator.dist(x2);  //aca evaluo en x2

			double deriv =  ( dist_x2/P2(x2)-dist_x/P2(x) )/ dx;    
			double deriv2 = ( dist_x2-dist_x )/ (P2(x)*dx) -  2.0*dist_x/P3(x);    
////////////////////////////

			double fSyn = P2(x)*Psyn*deriv;  // funcion a integrar

			L1   = L1 + fSyn*dx;

			//incremento(i,0,n-1,J1,J2,J3,L1);  //if (L1 > 0.0) 

			x = x*x_int ;

		}
	}

	double integral = L1;// (J1 + 2.0*J2 + 4.0*J3) / 3.0;

	double absorptionCoefficient = - P3(planck)*cLight2*integral/(8*pi*P2(E));

	double opacity = absorptionCoefficient*radius;

	return opacity;
}
