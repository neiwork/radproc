#include "neutronDecay.h"


#include <fmath/RungeKutta.h>
#include <fmath/physics.h>


double f_NeuDec(double Epi, double E, double mass, const Particle& neutron, const SpaceCoord& distCoord)         //funcion a integrar variable Epi
{
	
	double distNeutron;
	if (Epi < neutron.emin() || Epi > neutron.emax()){
		distNeutron = 0.0;
	}
	else{
		distNeutron = neutron.distribution.interpolate({ { 0, Epi } }, &distCoord); 
	}
	//double distNeutron = interpol(Epi,Ecreator,Ncreator,Ncreator.size()-1);

	double r = P2(mass/neutronMass);

	double x = E/Epi;
	double decayTime = neutronMeanLife*Epi/(neutronMass*cLight2);

	double Q = distNeutron/(Epi*decayTime); //*r*(1-x)/(Epi*x*P2(1-r)*decayTime);   

	return Q;   
}


double injNeutronDecay(double E, Particle& p, Particle& neutron, const SpaceCoord& distCoord)  
{

	double Emax = neutron.emax();  
	
	double mass = p.mass;
	
	//.Ecreator = neutron.eDim()->values; //el creator es neutron

	double injection = 0.0;

	double rpi = P2(mass/neutronMass);

	double sup = std::min(Emax,E/rpi);  // transformo la condicion de la heaviside en un limite superior

	injection = RungeKuttaSimple(E, sup, 
				[E,mass,&neutron,&distCoord](double Epi){
						return f_NeuDec(Epi,E,mass,neutron,distCoord);
					});
		
	return injection;

}

double f3(double gamma_3, Particle& n, const SpaceCoord& distCoord)
{
	double En = gamma_3 * neutronMass * cLight2;
	double Nn = (En > n.emin() && En < n.emax()) ? 
				n.distribution.interpolate({{0,En}},&distCoord) : 0.0;
	return (gamma_3 > 1.001) ? Nn / sqrt(gamma_3*gamma_3-1.0) / gamma_3 : 0.0;
}

double injElectronNeutronDecay(double Ee, Particle& n, const SpaceCoord& distCoord)
{
	double gamma_e = Ee/(electronMass*cLight2);
	double C = 0.5*(neutronMass/electronMass)/neutronMeanLife;
	double sup = 2.53;
	double inf = 1.001;
	size_t N1 = 20;
	double dgamma_1 = (sup-inf)/N1;
	double gamma_1 = inf;
	double sum = 0.0;
	for (size_t j=0;j<N1;j++) {
		double gamma_2_minus = gamma_e*gamma_1 - sqrt((gamma_1*gamma_1-1.0)*(gamma_e*gamma_e-1.0));
		double gamma_2_plus = gamma_e*gamma_1 + sqrt((gamma_1*gamma_1-1.0)*(gamma_e*gamma_e-1.0));
		double secondInt = RungeKuttaSimple(gamma_2_minus,gamma_2_plus,[&n,&distCoord]
							(double gamma_3) {return f3(gamma_3,n,distCoord);});
		sum += 0.614*gamma_1*P2(2.53-gamma_1)*dgamma_1*secondInt;
		gamma_1 += dgamma_1;
	}
	return C*sum;
}



/* la de arriba era protonNeutron, lo de abajo es de electrones, son lo mismo, cambia la masa

double fQelectron(double Epi, double E, const Particle& p, const SpaceCoord& distCoord)         //funcion a integrar variable Epi
{
	double mass = p.mass;
	
	double distNeutron = interpol(Epi,Ecreator,Ncreator,Ncreator.size()-1);

	double r = P2(mass/neutronMass);

	double x = E/Epi;
	double decayTime = neutronMeanLife*Epi/(neutronMass*cLight2);

	double Q = distNeutron/(Epi*decayTime); //*r*(1-x)/(Epi*x*P2(1-r)*decayTime);   

	return Q;   
}


double electronNeutron(double E, Vector Ncreator, Particle& particle, Particle& neutron)  
{
	ParticleType particleName = particle.type; 

	double Emax = 1.6e-12*pow(10.0,neutron.logEmax);  

	double mass = p.mass;

	DataInjection data;

	data.E        = E;
	data.mass     = particle.mass;
	data.Ncreator = Ncreator;
	data.Ecreator = neutron.energyPoints;

	double injection = 0.0;

	double rpi = P2(electronMass/neutronMass);

	double sup = std::min(Emax,E/rpi);  // transformo la condicion de la heaviside en un limite superior

	injection = RungeKuttaSimple(E, sup, fQelectron, &data);
		
	return injection;

}*/