#include "pairInjectionExact.h"


#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
//#include <algorithm>


double cAnnihilationExact(double e1, double E)
{

	double gamma = E/(electronMass*cLight2);
	double uno  = 1.0/e1;      //esta es la condicion sobre las energias de los
	                           //fotones para que puedan crear pares 
	double dos  = gamma+1.0-e1;
	return std::max(uno,dos);
}

double dAnnihilationExact(double e1, double targetPhotonEmax)
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

	if (c_mas == 0.0)
	{	
		H_mas = (P3(ecm)/12.0-ecm*d_mas/8.0)/pow(e1*e2, 1.5)+
			    (P3(ecm)/6.0+ecm/2.0+1.0/(4*ecm))/sqrt(e1*e2);
	}
	else 
	{
		if (c_mas > 0.0)
		{
			I_mas = pow(c_mas,-0.5)*log(ecm*sqrt(c_mas)+sqrt(e1*e2+c_mas*P2(ecm)));
		}
		else // (c_mas < 0)
		{
			I_mas = pow(-c_mas,-0.5)*asin(ecm*sqrt(-c_mas/(e1*e2)));
		}

		H_mas = -ecm*(d_mas/(e1*e2)+2.0/c_mas)/(8.0*sqrt(e1*e2+c_mas*P2(ecm)))+
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
		if (c_menos > 0.0)
		{
			I_menos = pow(c_menos,-0.5)*log(ecm*sqrt(c_menos)+sqrt(e1*e2+c_menos*P2(ecm)));
		}
		else  //(c_menos < 0)
		{
			I_menos = pow(-c_menos,-0.5)*asin(ecm*sqrt(-c_menos/(e1*e2)));
		}

		H_menos = -ecm*(d_menos/(e1*e2)+2.0/c_menos)/(8.0*sqrt(e1*e2+c_menos*P2(ecm)))+
			   0.25*(2.0-(e1*e2-1.0)/c_menos)*I_menos+
			   0.25*sqrt(e1*e2+c_menos*P2(ecm))*(ecm/c_menos+1.0/(ecm*e1*e2));
	}		

	return 0.25*sqrt(P2(Ep)-4.0*P2(ecm))+H_mas+H_menos;

}


double fAnnihilationExact(double e1, double e2, double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax) 
{ 
	//e1=Ega; e2=Eph
	double Erest = electronMass*cLight2;
	
	double gamma = E/Erest;

	double Ep = e1+e2;

	double photonDist_x(0.0), photonDist_y(0.0);
	double x = e1*Erest;
	double y = e2*Erest;
	if (x >= tpEmin && x <= tpEmax){
		photonDist_x = ntPh.interpolate({ { 0, x } }, &distCoord); 
	}
	double nph;
	if (y >= tpEmin && y <= tpEmax){
		photonDist_y = tpf.interpolate({ { 0, y } }, &distCoord); 
	}
	
	//double photonDist_x = tpf(e1*Erest);  //aca los evaluo en las Eph no normalizadas
	//double photonDist_y = tpf(e2*Erest);

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



//double pairInjectionExact(double E, Particle& particle, Particle& photon, fun1 tpf)
double pairInjectionExact(double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)
{	
	using std::bind; using namespace std::placeholders; // para _1, _2, etc.

	double Erest = electronMass*cLight2;
	
	double cte = 3.0*cLight*thomson*Erest/4.0;  // *Erest^2: de las dos dist de fotones
	                                              // /Erest: de la inyeccion de electrones

//normalizo todo a me*c^2

	double inf = tpEmin/Erest; 

	double sup = tpEmax/Erest;  

	double integral  = RungeKutta(inf,sup, 
		[E](double u) {return cAnnihilationExact(u, E); }, //max(tpEmin,cAnnihilation(u, E)); },  //limite inferior
		[tpEmax](double u) {return dAnnihilationExact(u, tpEmax); },							  //limite superior
		[E, tpf, ntPh, &distCoord, tpEmin, tpEmax](double u, double t) 
		{return fAnnihilationExact(u, t, E, ntPh, tpf, distCoord, tpEmin, tpEmax); });
		

	double emissivityA = cte*integral;

	return emissivityA;

}





/*

namespace {
	void incremento(int entero, int left, int right, double &A1, double &A2, double &A3, double R)
	{
			double floaT = entero;
			double aux = fmod(floaT,2.0);
			
			if(entero==left || entero==right) //si el entero coincide con los extremos
				{A1 += R ;}
			else if(aux==0)				//si el entero es par
				{A2 += R ;}
			else 
				{A3 += R ;}    //si el entero es impar

	}
}
//******************************************************************************************************
//******************************************************************************************************
//Emisividad de pares electrón-positrón por aniquilación fotón-fotón. Las
//fórmulas son del paper de Böttcher & Schlickeiser (1997), en la aproximación de
//Aharonian et al. (1983).

//double gamma_gamma_pair_emissivity(const double electronEnergy, 
//								   const double z, 
//								   const Particle_t &targetPhoton1,
//								   const Particle_t &targetPhoton2)
double pairInjection2(double E, 
                     Particle& particle, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)
{
	unsigned int energyIndex1;
	unsigned int energyIndex2;
	const unsigned int nPoints = 200;

	double electronEnergy = E;

	double targetPhoton2Emin = 1.0*1.6e-12;//pow(10.0,photonLogEmin)*1.6e-12;  //0.15e3*1.6e-12;
	double targetPhoton2Emax = 200.0e3*1.6e-12;//pow(10.0,photonLogEmax)*1.6e-12;  //150.0e3*1.6e-12;

	//double targetPhoton2Emin = pow(10.0,photonLogEmin)*1.6e-12;  //0.15e3*1.6e-12;
	//double targetPhoton2Emax = pow(10.0,photonLogEmax)*1.6e-12;  //150.0e3*1.6e-12;
	
	static Vector term(7, 0.0);

	double energy1, denergy1, energy1_int;;
	double epsilon1;
	double energy2, denergy2, energy2_int;; 
	double epsilon2;
	double epsilon12;

	double lowerLimit2;
	double lowerLimit1;
	const double upperLimit1 = 0.99*tpEmax;   
	const double upperLimit2 = 0.99*tpEmax;//targetPhoton2Emax; // energia max de los fotones de la corona                 
	const double gamma       = electronEnergy/(electronMass*cLight2);
	
	lowerLimit1 = std::max(1.01*electronEnergy, tpEmin);

	if(lowerLimit1 > upperLimit1)
	{
		return 0.0;
	} 

			//energyIntegral_1 = create_integral(nPoints, order_two);
			//energyIntegral_2 = create_integral(nPoints, order_two);

	term[0] = (3.0/32.0)*cLight*thomson/(electronMass*cLight2);

	double J1(0.0), J2(0.0), J3(0.0);

	energy1_int = pow((upperLimit1/lowerLimit1),(1.0/nPoints)) ;
			
	energy1 = lowerLimit1;
	
	for (energyIndex1 = 0; energyIndex1 < nPoints; ++energyIndex1)
	{
		denergy1  = energy1*(energy1_int-1.0);
		epsilon1  = energy1/(electronMass*cLight2);

				//energyIntegral_1.x[energyIndex1] = energy1;

				//double photonDist_x = data->tpf(x);
				//double photonDist_y = data->tpf(y);
	
		term[1] = tpf(energy1)/pow(epsilon1, 3.0);		
				//interpolate_double_data(energy1, z, targetPhoton1.distribution,targetPhoton1.energyPoints, targetPhoton1.zPoints,targetPhoton1.nEnergyPoints,targetPhoton1.nZPoints)/pow(epsilon1, 3.0);		

		lowerLimit2 = std::max(1.01*electronMass*cLight2*epsilon1/
									 (4.0*gamma*(epsilon1 - gamma)), 
									  targetPhoton2Emin);
		
		
		double K1(0.0), K2(0.0), K3(0.0);

		if(lowerLimit2 < upperLimit2)
		{
		
			energy2_int = pow((upperLimit2/lowerLimit2),(1.0/nPoints)) ;	

			energy2 = lowerLimit2;
			
			for(energyIndex2 = 0; energyIndex2 < nPoints; ++energyIndex2)
			{
				denergy2  = energy2*(energy2_int-1.0);
				
						//energy2  = logspace(energyIndex2, lowerLimit2, upperLimit2, nPoints);				
				epsilon2 = energy2/(electronMass*cLight2);

						//energyIntegral_2.x[energyIndex2] = energy2;

				epsilon12 = epsilon1*epsilon2;

				term[2] = tpf(energy2)/pow(epsilon2, 2.0);
							//interpolate_double_data(energy2, z, targetPhoton2.distribution,targetPhoton2.energyPoints, targetPhoton2.zPoints, targetPhoton2.nEnergyPoints,targetPhoton2.nZPoints)/pow(epsilon2, 2.0);	

			    term[3] = 4.0*P2(epsilon1)/(gamma*(epsilon1 - gamma))*
						  log(4.0*epsilon2*gamma*(epsilon1 - gamma)/epsilon1);
				
				term[4] = 8.0*epsilon12;

				term[5] = 2.0*(2.0*epsilon12 - 1.0)*
						  P2(epsilon1)/(gamma*(epsilon1 - gamma));

				term[6] = (1.0-1.0/(epsilon12))*P2(epsilon1)*P2(epsilon1)/
						   P2(gamma*(epsilon1 - gamma));

				double L1   = term[2]*(term[3] - term[4] + term[5] - term[6])*denergy2;
				
				if (L1 > 0.0) {incremento(energyIndex2,0,nPoints-1,K1,K2,K3,L1);}
				
							//energyIntegral_2.integrand[energyIndex2] = 
							//term[2]*(term[3] - term[4] + term[5] - term[6]);
				
				energy2 = energy2*energy2_int;				
					
			}  // for energyIndex2 		 
			
			double L2  = (K1+2.0*K2+4.0*K3)*denergy1*term[1]/3.0;

			if (L2 > 0.0) {incremento(energyIndex1,0,nPoints-1,J1,J2,J3,L2);}
				
							//sum[0] = integrate(energyIntegral_2);
							//energyIntegral_1.integrand[energyIndex1] = sum[0]*term[1];
		
		}
		
		energy1 = energy1*energy1_int;
		
	} // for energyIndex1 

	double result = (J1+2.0*J2+4.0*J3)*term[0]/3.0;
	
	return result;
	
				//sum[1] = integrate(energyIntegral_1);
				//return (term[0]*sum[1]);
				
	
}*/