#include "lossesPhotoHadronic.h"

#include "crossSectionInel.h"
#include <fmath\RungeKutta.h>
#include <fparameters\parameters.h>
#include <fmath\physics.h>
#include <algorithm>

double cPionPH(double u)   //limite inferior
{
	return pionThresholdPH; //145 MeV  threshold energy for pion production
}

double cPairPH(double u)   //limite inferior
{
	return pairThresholdPH;
}

double dPH(double u, double E, double mass)   //limite superior     
{

	return 2 * u*E / (mass*cLight2);   //d1 = 2*u*ep/(masa*cluz**2)
}

/*
//double fPHPion(double u,double t, double E, double mass, fun1 tpf)   //funcion a integrar
double fPHPion(double u, double t, double E, double mass, const ParamSpaceValues& tpf, const SpaceCoord& psc)   //funcion a integrar
{
	double Nph = tpf.interpolate({ { 0, u } }, &psc);
	//return tpf(u)*crossSectionPHPion(t)*inelasticityPHPion(t)*t/P2(u);
	return Nph*crossSectionPHPion(t)*inelasticityPHPion(t)*t / P2(u);
}

double fPHPair(double u, double t, double E, double mass, const ParamSpaceValues& tpf, const SpaceCoord& psc)   //funcion a integrar
{
	double Nph = tpf.interpolate({ { 0, u } }, &psc);
	return Nph
		*crossSectionBetheHeitler(t)*inelasticityBetheHeitler(t)*t / P2(u);
}

double lossesPhotoHadronic(double E, Particle& particle, const ParamSpaceValues& tpf, const SpaceCoord& psc, double phEmin, double phEmax)
{  //E=Ep

	double mass = particle.mass;

	double cte = 0.5*P2(mass*cLight2)*cLight;

	double b = phEmax;   //energia maxima de los fotones en erg

	double a1 = mass*cLight2*pionThresholdPH / (2 * E);
	double a2 = mass*cLight2*pairThresholdPH / (2 * E);

	double a_pi = std::max(a1, phEmin); //(a1,targetPhotonEmin); 
	double a_pa = std::max(a2, phEmin); //targetPhotonEmin);

	double integral = RungeKutta(a_pi, b, &cPionPH,
		[E, mass](double u) {
		return dPH(u, E, mass);
	}, [&](double u, double t) {
		return fPHPion(u, t, E, mass, tpf, psc);
	});


	integral += RungeKutta(a_pa, b, &cPairPH,
		[E, mass](double u) {
		return dPH(u, E, mass);
	}, [E, mass, tpf, &psc](double u, double t) {
		return fPHPair(u, t, E, mass, tpf, psc);
	});

	return cte*integral / E;
}
*/





double fPHPion(double u,double t, double E, double mass, fun1 tpf)   //funcion a integrar
{
	return tpf(u)*crossSectionPHPion(t)*inelasticityPHPion(t)*t/P2(u);
}

double fPHPair(double u,double t, double E, double mass, fun1 tpf)   //funcion a integrar
{
	double pepe = tpf(u);
	return (pepe)
		   *crossSectionBetheHeitler(t)*inelasticityBetheHeitler(t)*t/P2(u);
}

double lossesPhotoHadronic(double E, Particle& particle, fun1 tpf, double phEmin, double phEmax)  //E=Ep
{

	double mass = particle.mass;

	double cte	=	0.5*P2(mass*cLight2)*cLight;

	double b   = 10.0*phEmax;   //energia maxima de los fotones en erg

	double a1  = mass*cLight2*pionThresholdPH/(2*E);
	double a2  = mass*cLight2*pairThresholdPH/(2*E);

	double a_pi = std::max(a1,phEmin); 
	double a_pa = std::max(a2,phEmin);

	double integral = RungeKutta(a_pi, b, &cPionPH, 
	[E,mass](double u){
		return dPH(u,E,mass); 
	}, [E, mass, tpf](double u, double t){
		return fPHPion(u,t,E,mass,tpf);
	});

	integral += RungeKutta(a_pa, b, &cPairPH,
		[E, mass](double u){
		return dPH(u, E, mass);
	}, [E, mass, tpf](double u, double t){
		return fPHPair(u, t, E, mass, tpf);
	});

	return cte*integral/E;
}


