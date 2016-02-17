#include "pgammaPairInj.h"

#include "dataInjection.h"
#include <fmath\RungeKutta.h>
#include <fparameters\parameters.h>
#include <flosses\crossSectionInel.h>
#include <flosses\lossesPhotoHadronic.h>
#include <fmath\interpolation.h>

double fOmegaPHPair(double u,double t, void* voiddata)   //funcion a integrar
{
	DataInjection* data = (DataInjection*)voiddata;
	const double E = data->E;
	double mass = data->mass;

	return (data->tpf->call(u))
		     *crossSectionBetheHeitler(t)*t/P2(u);
}

double f_t_PHPair(double u,double t, void* voiddata)   //funcion a integrar
{
	DataInjection* data = (DataInjection*)voiddata;
	const double E = data->E;
	double mass = data->mass;

	double pepe = data->tpf->call(u);
	return (pepe)
		     *crossSectionBetheHeitler(t)*inelasticityBetheHeitler(t)*t/P2(u);
}



double omega_pair_PH(double E, Particle& particle, Fun* tpf)  //E=Ep
{
	DataInjection data;
	data.E = E;
	data.mass = particle.mass;
	data.tpf = tpf;

	double mass = particle.mass;
	double cte	= 0.5*P2(mass*cLight2)*cLight;

	double b   = 10.0*targetPhotonEmax;   //energia maxima de los fotones en erg

	double a1  = mass*P2(cLight)*pairThresholdPH/(2*E);

	double a = std::max(a1,targetPhotonEmin);

	double integral = RungeKutta(a,b,&cPairPH,&dPH,&fOmegaPHPair,&data);

	return cte*integral/P2(E);
}	

double t_pair_PH(double E, Particle& particle, Fun* tpf)  //E=Ep
{
	DataInjection data;
	data.E = E;
	data.mass = particle.mass;
	data.tpf  = tpf;

	double mass = particle.mass;
	double cte	=	0.5*P2(mass*cLight2)*cLight;

	double b   = 10*targetPhotonEmax;   //energia maxima de los fotones en erg

	double a1  = mass*P2(cLight)*pairThresholdPH/(2*E);

	double a = std::max(a1,targetPhotonEmin);

	double integral = RungeKutta(a,b,&cPairPH,&dPH,&f_t_PHPair,&data);

	return cte*integral/P2(E);
}


double pgammaPairInj2(double E, Vector Nproton, Particle& particle, Particle& proton, Fun* tpf)  
{
	double veinteE = 20.0*E;
	double protonDist = interpol(veinteE,proton.energyPoints,Nproton,Nproton.size()-1);
	//double protonDist = 3;

	double t_1   = t_pair_PH(veinteE, proton, tpf);     //esto no es lossesPH porque son perdidas solo del canal de produccion de piones
	double omega = omega_pair_PH(veinteE, proton, tpf);
	
	double emissivity;

	if (omega > 0.0)	{  //t_1 > 0.0 && 
		double averageInel = t_1/omega;

		double k1 = 0.2;

		double k2 = 0.6;

		double p1 = (k2-averageInel)/(k2-k1);

//		double pseda = 0.5;

		emissivity = 20.0*(2.0-1.5*p1)*omega*protonDist;
	}
	else	{
		emissivity = 0.0;
	}
	
	return emissivity;
}