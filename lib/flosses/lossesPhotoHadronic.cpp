#include "lossesPhotoHadronic.h"

#include "crossSectionInel.h"
#include "dataLosses.h"
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
	//DataLosses* data = (DataLosses*)voiddata;
	//const double E = data->E;
	//const double mass = data->mass;

	return 2*u*E/(mass*cLight2);   //d1 = 2*u*ep/(masa*cluz**2)
}

double fPHPion(double u,double t, double E, double mass, fun1 tpf)   //funcion a integrar
{
	//DataLosses* data = (DataLosses*)voiddata;
	//const double E = data->E;
	//const double mass = data->mass;

	return tpf(u)*crossSectionPHPion(t)*inelasticityPHPion(t)*t/P2(u);
}

double fPHPair(double u,double t, double E, double mass, fun1 tpf)   //funcion a integrar
{
	//DataLosses* data = (DataLosses*)voiddata;
	//const double E = data->E;
	//const double mass = data->mass;

	double pepe = tpf(u);
	return (pepe)
		   *crossSectionBetheHeitler(t)*inelasticityBetheHeitler(t)*t/P2(u);
}

double lossesPhotoHadronic(double E, Particle& particle, fun1 tpf)  //E=Ep
{

	double mass = particle.mass;

	double cte	=	0.5*P2(mass*cLight2)*cLight;

	double b   = 10.0*targetPhotonEmax;   //energia maxima de los fotones en erg

	double a1  = mass*cLight2*pionThresholdPH/(2*E);
	double a2  = mass*cLight2*pairThresholdPH/(2*E);

	double a_pi = std::max(a1,targetPhotonEmin); 
	double a_pa = std::max(a2,targetPhotonEmin);

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



//double cPH(double u, void*)   //limite inferior
//{
//	return pionThresholdPH; //145 MeV  threshold energy for pion production
//}
//
//double dPH(double u, void* voiddata)   //limite superior     
//{
//	DataLosses* data = (DataLosses*)voiddata;
//	const double E = data->E;
//	const double mass = data->mass;
//
//	return 2*u*E/(mass*cLight2);   //d1 = 2*u*ep/(masa*cluz**2)
//}
//
//double fPHPion(double u,double t, void* voiddata)   //funcion a integrar
//{
//	DataLosses* data = (DataLosses*)voiddata;
//	const double E = data->E;
//	const double mass = data->mass;
//
//	return (data->tpf(u))
//		     *crossSectionPHPion(t)*inelasticityPHPion(t)*t/P2(u);
//}
//
//double fPHPair(double u,double t, void* voiddata)   //funcion a integrar
//{
//	DataLosses* data = (DataLosses*)voiddata;
//	const double E = data->E;
//	const double mass = data->mass;
//	fun1 targetPhotonField = data->tpf;
//
//	double Erep = mass*cLight2;
//	double r    = t*P2(Erep)/(4*u*E*(E-t));
//	double pepe = targetPhotonField(u);
//
//	return (pepe) //data->tpf(u))
//		   *crossSectionBetheHeitler(t)*inelasticityBetheHeitler(t)*t/P2(u);
//}
//
//double lossesPhotoHadronic(double E, void* data, Particle& particle)  //E=Ep
//{
//
//	double mass = particle.mass;
//
//	double cte	=	0.5*P2(mass*cLight2)*cLight;
//
//	double b   = 10*targetPhotonEmax;   //energia maxima de los fotones en erg
//
//	double a1  = mass*P2(cLight)*pionThresholdPH/(2*E);
//
//	double a = std::max(a1,targetPhotonEmin);
//
//	double integral = //RungeKutta(a,b,&cPH,&dPH,&fPHPion,data) ;
//		              + RungeKutta(a,b,&cPH,&dPH,&fPHPair,data);
//
//	return cte*integral/E;
//}