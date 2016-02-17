#include "pairInjection.h"

#include "DataInjection.h"
#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
//#include <fmath\interpolation.h>
#include <fmath\physics.h>
#include <algorithm>

double cAnnihilation(double x,double E, double mass)
{
	//DataInjection* data = (DataInjection*)voiddata;
	//double E = data->E;      
	//double electronMass = data->mass;

	const double Erep = mass*cLight2;
	
	double inf = x*P2(Erep)/(4*(E)*(x-(E)));
//targetPhotonEmin
	double result = std::max(inf,0.2e3*1.6e-12);     //elijo el maximo entre el limite inferior y la minima energía de
	                                 //0.15 o 0.2                   //los fotones blanco
	return result;   //GUARDA que aca iba result
}

double dAnnihilation(double x )
{
	return targetPhotonEmax;                //este es el infinito del limite superior para la integral en Eph
}

double fAnnihilation(double x, double y, double E, fun1 tpf)   //x=Ega; y=Eph
{ 
	//DataInjection* data = (DataInjection*)voiddata;
	//double E = data->E;     
	//double electronMass = data->mass;
//	Vector& NgaP = data->Ncreator;
//	Vector& EgaP = data->Ecreator;

	double photonDist_x = tpf(x);
	double photonDist_y = tpf(y);

//	double lum = interpol(x,EgaP,NgaP,NgaP.size()-1);

//	double dist = lum*radius/(volume*P2(x)*cLight);  //lo cambie porque estaba mal  asi o lum/(4*pi*P2(R)*cLight*P2(E))  ???


	const double Erest = electronMass*cLight2;
                       //VER!!!!!!!!!!!!!!
	double result  =   ((photonDist_x*photonDist_y)/(P2(y)*P3(x)))*       //+nphBlackBody  lo saco porque Ega > 10^15 eV para crear
		                                                                  //pares a partir de la interaccion con los fotones del disco
					   ((4*P2(x)/(E*(x-E)))*log(4*E*y*(x-E)/(P2(Erest)*x))
					   -8*x*y/(P2(Erest))+2*(2*x*y-P2(Erest))*P2(x/Erest)/(E*(x-E))
					   -(1-P2(Erest)/(x*y))*pow(x,4)/(P2(E*(x-E))));
//&& (y > 2*P2(Erest)/x)
	return ((y<Erest) && (y > 2*P2(Erest)/x)  ) ? result : 0.0;   //pido que de algo solo si epsilon < Erep   //&& Erest<x
}


double pairInjectionFlor(double E, Particle& particle, Particle& photon, fun1 tpf)
{	
	using std::bind; using namespace std::placeholders; // para _1, _2, etc.

	double cte = 3.0*cLight*thomson*pow((particle.mass*cLight2),4)/32.0;

	double sup = 1.0; // 2e9*1.6e-12;  //este es el infinito del limite superior para la integral en Egamma
                                         //como Eph_min =0.15 keV --> Ega_max = 2 GeV 														

	double inf = E;  //Ega_min < Ee_min  --> la condicion esta asegurada
	
	//DataInjection data;

	//data.E        = E;
	//data.mass     = particle.mass;
	//data.tpf      = tpf;

	double integral  = RungeKutta(inf,sup,
		bind(cAnnihilation,_1,E,particle.mass),
		dAnnihilation,
		bind(fAnnihilation,_1,_2,E,tpf));

	double emissivityA = cte*integral;

	return emissivityA;

}


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
double pairInjection(double E, 
                     Vector Nphoton, Particle& particle, Particle& photon, fun1 tpf)
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
	const double upperLimit1 = 0.99*targetPhotonEmax;   
	const double upperLimit2 = targetPhoton2Emax; // energia max de los fotones de la corona                 
	const double gamma       = electronEnergy/(electronMass*cLight2);
	
	lowerLimit1 = std::max(1.01*electronEnergy, targetPhotonEmin);

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
				
	
}

