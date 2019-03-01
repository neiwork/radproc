#include "pgammaPionInj.h"


#include <fmath\RungeKutta.h>
#include <fparameters\parameters.h>
#include <flosses\crossSectionInel.h>
#include <flosses\lossesPhotoHadronic.h>
#include <fmath\interpolation.h>
#include <algorithm>

double fOmegaPHPion(double u,double t, double E, double mass, fun1 tpf)   //funcion a integrar
{

	return (tpf(u))
		     *crossSectionPHPion(t)*t/P2(u);
}

double f_t_PHPion(double u,double t, fun1 tpf)   //funcion a integrar
{

	double pepe = tpf(u);
	return (pepe)
		     *crossSectionPHPion(t)*inelasticityPHPion(t)*t/P2(u);
}


double omegaPH(double E, const Particle& particle, fun1 tpf, double tpEmin, double tpEmax)  //E=Ep
{
	using std::bind; using namespace std::placeholders; // para _1, _2, etc.

	double mass = particle.mass;
	double cte	= 0.5*P2(mass*cLight2)*cLight;

	double b   = 10.0*tpEmax;   //energia maxima de los fotones en erg

	double a1  = mass*P2(cLight)*pionThresholdPH/(2*E);

	double a = std::max(a1,tpEmin);

	double integral = RungeKutta(a, b, &cPionPH, 
		bind(dPH,_1, E, mass),
		bind(fOmegaPHPion,_1,_2,E,mass,tpf)
	);

	return cte*integral/P2(E);
}	

double t_pion_PH(double E, const Particle& particle, fun1 tpf, double tpEmin, double tpEmax)  //E=Ep
{

	double mass = particle.mass;
	double cte	=	0.5*P2(mass*cLight2)*cLight;

	double b   = 10*tpEmax;   //energia maxima de los fotones en erg

	double a1  = mass*P2(cLight)*pionThresholdPH/(2*E);

	double a = std::max(a1,tpEmin);

	double integral = RungeKutta(a, b, &cPionPH, [E, mass](double u){return dPH(u, E, mass); }, [tpf](double u, double t){
		return f_t_PHPion(u,t,tpf); 
	});

	return cte*integral/P2(E);
}


//double pgammaPionInj(double E, Vector Nproton, Particle& particle, Particle& proton, fun1 tpf)  
double pgammaPionInj(double E, const Particle& creator,	const SpaceCoord& psc, fun1 tpf, double tpEmin, double tpEmax)
{
	double cincoE = 5.0*E;
	double protonDist = creator.distribution.interpolate({ { 0, cincoE } }, &psc); //proton.dist(cincoE);

	double t_1   = t_pion_PH(cincoE, creator, tpf, tpEmin, tpEmax);     //esto no es lossesPH porque son perdidas solo del canal de produccion de piones
	double omega = omegaPH(cincoE, creator, tpf, tpEmin, tpEmax);
	
	double emissivity;

	if (t_1 > 0.0 && omega > 0.0)	{
		double averageInel = t_1/omega;

		double k1 = 0.2;

		double k2 = 0.6;

		double p1 = (k2-averageInel)/(k2-k1);
	
		double nChargedPion = 2.0-1.5*p1;

		emissivity = 5.0*nChargedPion*omega*protonDist;
	}
	else	{
		emissivity = 0;
	}
	
	return emissivity;
}


//no funciono este, me dan muchos ceros a baja energia

//double pgammaPairInj(double E, Vector Nproton, Particle& particle, Particle& proton, fun1 tpf)  
//{
//	double veinteE = 20.0*E;
//	double protonDist = interpol(veinteE,proton.energyPoints,Nproton,Nproton.size()-1);
//	//double protonDist = 3;
//
//	double t_1   = t_pion_PH(veinteE, proton, tpf);     //esto no es lossesPH porque son perdidas solo del canal de produccion de piones
//	double omega = omegaPH(veinteE, proton, tpf);
//
//	double emissivity;
//
//	if (omega > 0.0)	{  //t_1 > 0.0 && 
//		double averageInel = t_1/omega;
//
//		double k1 = 0.2;
//
//		double k2 = 0.6;
//
//		double p1 = (k2-averageInel)/(k2-k1);
//
////		double pseda = 0.5;
//
//		emissivity = 20.0*(2.0-1.5*p1)*omega*protonDist;
//	}
//	else	{
//		emissivity = 0.0;
//	}
//	
//	return emissivity;
//}