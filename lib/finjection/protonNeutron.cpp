#include "protonNeutron.h"

#include "dataInjection.h"
#include <fmath\RungeKutta.h>
#include <fmath\physics.h>
#include <fmath\interpolation.h>
#include <algorithm>

//
//double fQproton(double Epi, void* voiddata)         //funcion a integrar variable Epi
//{
//	DataInjection* data = (DataInjection*)voiddata;
//	double E = data->E;      //E=Emu
//	double mass   = data->mass;
//	Vector& Ncreator = data->Ncreator;
//	Vector& Ecreator = data->Ecreator;
//
//	double distNeutron = interpol(Epi,Ecreator,Ncreator,Ncreator.size()-1);
//
//	double r = P2(mass/neutronMass);
//
//	double x = E/Epi;
//
//	double decayTime = neutronMeanLife*Epi/(neutronMass*cLight2);
//
//	double Q = distNeutron/(Epi*decayTime); //*r*(1-x)/(Epi*x*P2(1-r)*decayTime);   
//
//	
//	return Q;   
//}
//
//
//
//double protonNeutron(double E, Vector Ncreator, Particle& particle, Particle& neutron)  
//{
//	ParticleType particleName = particle.type; 
//
//	double Emax = 1.6e-12*pow(10.0,neutron.logEmax);  
//
//	DataInjection data;
//
//	data.E        = E;
//	data.mass     = particle.mass;
//	data.Ncreator = Ncreator;
//	data.Ecreator = neutron.eDim()->values;
//
//	double injection = 0.0;
//
//	double rpi = P2(protonMass/neutronMass);
//
//	double sup = std::min(Emax,E/rpi);  // transformo la condicion de la heaviside en un limite superior
//
//	injection = RungeKuttaSimple(E, sup, fQproton, &data);
//		
//	return injection;
//
//}
