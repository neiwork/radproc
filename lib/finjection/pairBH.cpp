#include "pairBH.h"


#include "dataInjection.h"
#include <fmath\RungeKutta.h>
#include <fparameters\parameters.h>
#include <flosses\lossesPhotoHadronic.h>
#include <flosses\crossSectionInel.h>
#include <fmath\interpolation.h>
#include <algorithm>


//double cBH(double u, void*)   //limite inferior
//{
//	return pairThresholdPH; //2*me*c^2  threshold energy for pair production
//} esto es igual a cPairPH

//double dBH(double u, void* voiddata)   //limite superior     
//{
//	DataInjection* data = (DataInjection*)voiddata;
//	const double E = data->E;
//	const double mass = data->mass;
//
//	return 2*u*E/(mass*cLight2);   //d1 = 2*u*ep/(masa*cluz**2)
//}  esto es igual a dPH

double dBH(double u, double E, double mass)   //limite superior     
{
	//DataInjection* data = (DataInjection*)voiddata;
	//const double E = data->E;
	//const double mass = data->mass;

	return 2*u*E/(mass*cLight2);   //d1 = 2*u*ep/(masa*cluz**2)
}


double fOmegaBH(double u,double t, double E, double mass, fun1 tpf)   //funcion a integrar
{
	//DataInjection* data = (DataInjection*)voiddata;
	//const double E = data->E;
	//double mass = data->mass;

	return (tpf(u))   //esto es el n(epsilon)
		     *crossSectionBetheHeitler(t)*t/P2(u);
}


double omegaBH(double E, Particle& particle, fun1 tpf)  //E=Ep
{
	//DataInjection data;
	//data.E = E;
	//data.mass = particle.mass;
	//data.tpf = tpf;

	double mass = particle.mass;
	double cte	=	0.5*P2(mass*cLight2)*cLight;

	double b   = 10*targetPhotonEmax;   //energia maxima de los fotones en erg

	double a1  = mass*cLight2*pairThresholdPH/(2*E);

	double a = std::max(a1,targetPhotonEmin);

	double integral = RungeKutta(a, b, &cPairPH, [E,mass](double u){
		return dPH(u, E, mass); 
	}, [E, mass, tpf](double u, double t){
		return fOmegaBH(u, t, E, mass, tpf);
	});

	return cte*integral/P2(E);
}	

double pairBH(double E, Particle& particle, Particle& proton, fun1 tpf)  
{
	double evalE = E*proton.mass/particle.mass;
	double protonDist = proton.dist(evalE);// interpol(evalE, proton.energyPoints, Nproton, Nproton.size() - 1);
	//double protonDist = 3;

	double omega = omegaBH(evalE, proton, tpf);  //de pares
	
	double emissivity = 2.0*proton.mass*omega*protonDist/particle.mass;
	
	return emissivity;
}



/*double fBH(double u,double t, void* voiddata)   //funcion a integrar
{
	DataInjection* data = (DataInjection*)voiddata;
	const double E = data->E;
	double mass = data->mass;
	Vector& Ncreator = data->Ncreator;
	Vector& Ecreator = data->Ecreator;

	double protonDist = interpol(u,Ecreator,Ncreator,Ncreator.size()-1);

	double pepe = data->tpf(t);
	return (pepe)
		     *crossSectionBetheHeitler(t)*pepe*protonDist; //*inelasticityBetheHeitler(t)*t/P2(u);
}

double pairBH(double E, Vector Nproton, Particle& particle, Particle& proton, fun1 tpf)  
{
	
	DataInjection data;
	data.E = E;
	data.mass = particle.mass;
	data.tpf  = tpf;
	data.Ecreator = proton.energyPoints;
	data.Ncreator = Nproton;

	double mass = particle.mass;
	double cte	= cLight/(mass*cLight2);

	double b = pow(10,proton.logEmax)*1.6e-12;   //energia maxima de los protones en erg

	double a = pow(10,proton.logEmin)*1.6e-12;   //energia minima de los protones en erg

	double integral = RungeKutta(a,b,&cBH,&dBH,&fBH,&data);

	return cte*integral;
}*/ 