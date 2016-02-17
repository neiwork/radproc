#include "pairInjectionExact.h"

#include "DataInjection.h"
#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\interpolation.h>
#include <fmath\physics.h>
#include <algorithm>


double cAnnihilationExact(double e1, double E, double mass)
{

	//DataInjection* data = (DataInjection*)voiddata;
	//double E = data->E;      
	//double electronMass = data->mass;  

	double gamma = E/(mass*cLight2);

	
	double uno  = 1.0/e1;      //esta es la condicion sobre las energias de los
	                           //fotones para que puedan crear pares 
	double dos  = gamma+1.0-e1;   

	double result1 = std::max(uno,dos);

	return result1;
}

double dAnnihilationExact(double e1 )
{
	double result = targetPhotonEmax/(electronMass*cLight2); //este es el infinito del limite superior 
	return result;                                       //para la integral en Eph  
}



double integrando (double ecm, double e1, double e2, double gamma)
{
	double Ep = e1+e2;
//////////definición de H_mas////////////////////////////////////
	double H_mas;
	double I_mas;
	
	double c_mas   = P2(e1-gamma)-1.0; 

	double d_mas   = P2(e1)+e1*e2+gamma*(e2-e1);

	if (c_mas == 0)
	{	
		H_mas = (P3(ecm)/12.0-ecm*d_mas/8.0)/pow(e1*e2, 1.5)+
			    (P3(ecm)/6.0+ecm/2.0+1.0/(4*ecm))/sqrt(e1*e2);
	}
	else 
	{
		if (c_mas > 0)
		{
			I_mas = pow(c_mas,-0.5)*log(ecm*sqrt(c_mas)+sqrt(e1*e2+c_mas*P2(ecm)));
		}
		else // (c_mas < 0)
		{
			I_mas = pow(-c_mas,-0.5)*asin(ecm*sqrt(-c_mas/(e1*e2)));
		}

		H_mas = -ecm*(d_mas/(e1*e2)+2.0/c_mas)/(8*sqrt(e1*e2+c_mas*P2(ecm)))+
			   0.25*(2.0-(e1*e2-1.0)/c_mas)*I_mas+
			   0.25*sqrt(e1*e2+c_mas*P2(ecm))*(ecm/c_mas+1.0/(ecm*e1*e2));
	}	
////////////////definición de H_menos/////////////////////////////////
	double H_menos;
	double I_menos;

	double c_menos = P2(e2-gamma)-1.0; 

	double d_menos = P2(e2)+e1*e2-gamma*(e2-e1);

	if (c_menos == 0)
	{	
		H_menos = (P3(ecm)/12.0-ecm*d_menos/8.0)/pow(e1*e2, 1.5)+
			      (P3(ecm)/6.0+ecm/2.0+1.0/(4*ecm))/sqrt(e1*e2);
	}
	else 
	{
		if (c_menos > 0)
		{
			I_menos = pow(c_menos,-0.5)*log(ecm*sqrt(c_menos)+sqrt(e1*e2+c_menos*P2(ecm)));
		}
		else  //(c_menos < 0)
		{
			I_menos = pow(-c_menos,-0.5)*asin(ecm*sqrt(-c_menos/(e1*e2)));
		}

		H_menos = -ecm*(d_menos/(e1*e2)+2.0/c_menos)/(8*sqrt(e1*e2+c_menos*P2(ecm)))+
			   0.25*(2.0-(e1*e2-1.0)/c_menos)*I_menos+
			   0.25*sqrt(e1*e2+c_menos*P2(ecm))*(ecm/c_menos+1.0/(ecm*e1*e2));
	}		

	return 0.25*sqrt(P2(Ep)-4*P2(ecm))+H_mas+H_menos;

}


double fAnnihilationExact(double e1, double e2, double E, double mass, fun1 tpf)   //e1=Ega; e2=Eph
{ 
	//DataInjection* data = (DataInjection*)voiddata;
	//double E = data->E;     
	//double electronMass = data->mass;

	double Erest = mass*cLight2;
	
	double gamma = E/Erest;

	double Ep = e1+e2;

	double photonDist_x = tpf(e1*Erest);  //aca los evaluo en las Eph no normalizadas
	double photonDist_y = tpf(e2*Erest);

	double ecm_a = sqrt(0.5*(gamma*(Ep-gamma)+1.0+sqrt(P2(gamma*(Ep-gamma)+1.0)-P2(Ep))));
	double ecm_d = sqrt(0.5*(gamma*(Ep-gamma)+1.0-sqrt(P2(gamma*(Ep-gamma)+1.0)-P2(Ep))));

	double ecmU = std::min(sqrt(e1*e2),ecm_a);
	double ecmL = std::max(1.0,ecm_d);

	double result  =   ((photonDist_x*photonDist_y)/P2(e1*e2))   
					   *(integrando(ecmU,e1,e2,gamma)-integrando(ecmL,e1,e2,gamma));
	
	//&& (e2 > 2*P2(Erest)/e1) esta condicion ya la agregue en el limite inferior sobre e2
	//return ((e2<Erest)  ) ? result : 0.0;
	
	return result;
}



double pairInjectionExact(double E, Particle& particle, Particle& photon, fun1 tpf)
{	
	using std::bind; using namespace std::placeholders; // para _1, _2, etc.

	double Erest = electronMass*cLight2;
	
	double cte = 3.0*cLight*thomson*Erest/4.0;  // *Erest^2: de las dos dist de fotones
	                                              // /Erest: de la inyeccion de electrones


	//normalizo todo a me*c^2

	double inf = targetPhotonEmin/Erest; 

	double sup = targetPhotonEmax/Erest;  											

	/*DataInjection data;

	data.E        = E;
	data.mass     = electronMass;
	data.tpf      = tpf;*/

	double integral  = RungeKutta(inf,sup,
		bind(cAnnihilationExact,_1,E,electronMass),
		dAnnihilationExact,
		bind(fAnnihilationExact,_1,_2,E,electronMass,tpf));

	double emissivityA = cte*integral;

	return emissivityA;

}