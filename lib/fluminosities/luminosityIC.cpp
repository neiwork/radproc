#include "luminosityIC.h"

#include <finjection\dataInjection.h>
#include <fmath\RungeKutta.h>
#include <fmath\interpolation.h>
#include <flosses\crossSectionInel.h>
#include <fparameters\parameters.h>
#include <fmath\physics.h>
#include <algorithm>


double cICemi(double u, double E, double mass)   //limite inferior
{
	//DataInjection* data = (DataInjection*)voiddata;
	//const double E = data->E;
	//double mass = data->mass;

	double Erep = mass*cLight2;

	double condition = E*u-P2(u);  //Ega*Ee-Ee^2

	double inf;

	if (condition > 0.0){
		inf = E*P2(Erep)/(4.0*condition);
	}
	else {
		inf = -E*P2(Erep)/(4.0*condition);
	}

	return std::max(targetPhotonEmin,inf);  //puse la condicion Ega < s*Ee/1+s
}

double dICemi(double u, double E)         //limite superior   
{
	//DataInjection* data = (DataInjection*)voiddata;
	//const double E = data->E;             

	return E;  //esta es la condicion epsilon < Ega	                                 
}

double fICemi(double u,double t, double E,const Particle& creator, fun1 tpf)   //funcion a integrar  u=Ee
{                                                                      //t=epsilon
	//DataInjection* data = (DataInjection*)voiddata;
	//const double E = data->E;                //E=Ega; L=L(Ega)
	//double mass = data->mass;
	//Vector& Ncreator = data->Ncreator;
	//Vector& Ecreator = data->Ecreator;

	double distCreator = creator.dist(u);// interpol(u, Ecreator, Ncreator, Ncreator.size() - 1);

	double Erep = creator.mass*cLight2;
	  

	double s = 4*u*t/P2(Erep);       //equivalente al gamma

	double r = E/(s*u*(1-E/u));   //equivalente al q

	double pepe = (tpf(t));

	double function = distCreator*(pepe/t)
            			*(2*r*log(r)+(1+2*r)*(1-r)+(1-r)*P2(r*s)/(2*(1+r*s)))/P2(u);
	
	double condition = s*u/(1+s);

	return function;

	//estas condiciones las transforme en dos limites de integracion
	//if (E<condition && t<E) return function;
	//else return 0.0;   //si se cumple las condiciones da funcion, sino cero
}

double luminosityIC(double E, const Particle& creator, fun1 tpf)
{
	//DataInjection data;

	//data.E = E;
	//data.mass = creator.mass;
	//data.Ncreator = Ncreator;
	//data.Ecreator = creator.energyPoints;
	//data.tpf      = tpf;
	double mass = creator.mass;
	
	double cte  = 3.0*crossSectionThomson(creator.mass)*P2(creator.mass)*pow(cLight,5)/4;

	double integral = RungeKutta(creator.emin(), creator.emax(),
		[E,mass](double u){
			return cICemi(u,E,mass);
		}, 
		[E](double u){
			return dICemi(u,E);
		}, 
		[E,&creator,tpf](double u, double t){
			return fICemi(u, t,E,creator,tpf); 
		});    //le asigno a la variable integral el resultado de la integracion   

	double luminosity = integral*cte*volume*P2(E);

	return luminosity;

	}

