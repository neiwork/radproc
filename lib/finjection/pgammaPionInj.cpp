#include "pgammaPionInj.h"


#include <fmath/RungeKutta.h>
#include <fparameters/parameters.h>
#include <flosses/crossSectionInel.h>
#include <flosses/lossesPhotoHadronic.h>
//#include <fmath/interpolation.h>
//#include <algorithm>

double fOmegaPHPion(double u,double t, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)//  fun1 tpf)   //funcion a integrar
{
	double nph;
	if (u < tpEmin || u> tpEmax){
		nph = 0.0;
	}
	else{
		nph = tpf.interpolate({ { 0, u } }, &distCoord); //tpf(u);
	}
	
	return nph*crossSectionPHPion(t)*t/P2(u);
}

double f_t_PHPion(double u,double t, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)//  fun1 tpf)   //funcion a integrar
{

	double nph;
	if (u < tpEmin || u> tpEmax){
		nph = 0.0;
	}
	else{
		nph = tpf.interpolate({ { 0, u } }, &distCoord); //tpf(u);
	}
	
	return nph*crossSectionPHPion(t)*inelasticityPHPion(t)*t/P2(u);
}


double omegaPH(double E, const Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& distCoord,
				double tpEmin, double tpEmax)  //E=Ep      fun1 -> const ParamSpaceValues& tpf
{
	//using std::bind; using namespace std::placeholders; // para _1, _2, etc.
	double g = E / (p.mass*cLight2);
	double cte	= 0.5*P2(p.mass*cLight2)*cLight;

	double b   = 10.0*tpEmax;   //energia maxima de los fotones en erg

	double a1  = p.mass*P2(cLight)*pionThresholdPH/(2*E);

	double a = std::max(a1,tpEmin);

	double integral = RungeKutta(a, b, &cPionPH, [E, &p](double u) {return dPH(u, E, p.mass); },  //[g](double u) {return dPH(u, g); },
		[tpf, &distCoord, tpEmin, tpEmax](double u, double t) {	return fOmegaPHPion(u, t, tpf, distCoord, tpEmin, tpEmax); });
		//bind(dPH,_1, E, mass),
		//bind(fOmegaPHPion,_1,_2,E,mass,tpf)

	return cte*integral/P2(E);
}

double integSigmaPG(double Eph, double g)
{
	double var = 2.0*Eph*g;
	if (var > 3.2e-4 && var < 8.0e-4)
		return 3.4e-28 * (var*var-1.03e-7);
	else if (var > 8.0e-4)
		return 1.83e-34 + 1.2e-28 * (var*var-6.42e-7);
	else
		return 0.0;
}

double integKSigmaPG(double Eph, double g)
{
	double var = 2.0*Eph*g;
	if (var > 3.2e-4 && var < 8.0e-4)
		return 0.2 * 3.4e-28 * (var*var-1.03e-7);
	else if (var > 8.0e-4)
		return 0.2 * 1.83e-34 + 0.6 * 1.2e-28 * (var*var-6.42e-7);
	else
		return 0.0;
}

double omegaPHsimple(double E, const Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& distCoord,
				double tpEmin, double tpEmax)
{
	double g = E / (p.mass*cLight2);
	double inf = pionThresholdPH/(2.0*g);
	return 0.25*cLight/(g*g) * integSimpsonLog(inf,tpEmax,[g,&tpf,&distCoord](double Eph)
				{
					double nPh = tpf.interpolate({{0,Eph}},&distCoord);
					return nPh/(Eph*Eph) * integSigmaPG(Eph,g);
				},30);
}	

double t_pion_PH(double E, const Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, 
					double tpEmin, double tpEmax)
{
	double g = E/(p.mass*cLight2);
	double cte	=	0.5*P2(p.mass*cLight2)*cLight;

	double b   = 10*tpEmax;   //energia maxima de los fotones en erg

	double a1  = p.mass*P2(cLight)*pionThresholdPH/(2*E);

	double a = std::max(a1,tpEmin);

 	double integral = RungeKutta(a, b, &cPionPH, [E, &p](double u) {return dPH(u, E, p.mass); },  //[g](double u) {return dPH(u, g); },
		[tpf, &distCoord, tpEmin, tpEmax](double u, double t){ return f_t_PHPion(u, t, tpf, distCoord, tpEmin,tpEmax); });

	return cte*integral/P2(E);
}

double t_pion_PHsimple(double E, const Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, 
					double tpEmin, double tpEmax)
{
	double g = E/(p.mass*cLight2);
	double inf = pionThresholdPH/(2.0*g);
	return 0.25*cLight/(g*g) * integSimpsonLog(inf,tpEmax,[g,&tpf,&distCoord](double Eph)
				{
					double nPh = tpf.interpolate({{0,Eph}},&distCoord);
					return nPh/(Eph*Eph) * integSigmaPG(Eph,g);
				},30);
	
}

//double pgammaPionInj(double E, Vector Nproton, Particle& particle, Particle& proton, fun1 tpf)  
double pgammaPionInj(double E, const Particle& creator,	const ParamSpaceValues& tpf, 
						const SpaceCoord& distCoord, double tpEmin, double tpEmax)
{
	double cincoE = 5.0*E;
	double protonDist = (cincoE > creator.emin() && cincoE < creator.emax()) ? 
					creator.distribution.interpolate({ { 0, cincoE } }, &distCoord) : 0.0; //proton.dist(cincoE);

	//double t_1   = t_pion_PHsimple(cincoE, creator, tpf, distCoord, tpEmin, tpEmax);     //esto no es lossesPH porque son perdidas solo del canal de produccion de piones
	//double omega = omegaPHsimple(cincoE, creator, tpf, distCoord, tpEmin, tpEmax);

	double t_1 = t_pion_PH(cincoE, creator, tpf, distCoord, tpEmin, tpEmax);     //esto no es lossesPH porque son perdidas solo del canal de produccion de piones
	double omega = omegaPH(cincoE, creator, tpf, distCoord, tpEmin, tpEmax);
	
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
