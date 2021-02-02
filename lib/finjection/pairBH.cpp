#include "pairBH.h"

//#include <gsl/gsl_math.h>
#include <fmath/RungeKutta.h>
#include <fparameters/parameters.h>
//#include <flosses/lossesPhotoHadronic.h>
#include <flosses/crossSectionInel.h>
//#include <algorithm>


double dBH(double u, double E, double mass)   //limite superior     
{
	return 2*u*E/(mass*cLight2);   //d1 = 2*u*ep/(masa*cluz**2)
}


double fOmegaBH(double u,double t, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)//  fun1 tpf)   //funcion a integrar
{
	double nph = 0.0;;
	if (u > tpEmin && u < tpEmax){
		nph = tpf.interpolate({ { 0, u } }, &distCoord); //tpf(u);
	}
	return nph*crossSectionBetheHeitler(t)*t/P2(u);
}



double omegaBH(double E, const Particle& particle, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)  //E=Ep      fun1 -> const ParamSpaceValues& tpf
{
	//using std::bind; using namespace std::placeholders; // para _1, _2, etc.

	double mass = particle.mass;
	double cte	= 0.5*P2(mass*cLight2)*cLight;

	double b   = 10.0*tpEmax;   //energia maxima de los fotones en erg

	double a1  = mass*cLight2*pairThresholdPH/(2.0*E);

	double a = std::max(a1,tpEmin);

	double integral = RungeKutta(a, b, 
					  [](double u) {return pairThresholdPH;},
					  [E, mass](double u) {return dBH(u, E, mass); },
		[tpf, &distCoord, tpEmin, tpEmax](double u, double t) {	return fOmegaBH(u, t, tpf, distCoord, tpEmin, tpEmax); });
		//bind(dPH,_1, E, mass),
		//bind(fOmegaPHPion,_1,_2,E,mass,tpf)

	return cte*integral/P2(E);
}

double omegaBH2(double e, const Particle& c, const ParamSpaceValues& tpf, const SpaceCoord& distCoord,
					double phEmin, double phEmax)
{
	double inf = pairThresholdPH*protonMass*cLight2/(2.0*e);
	double integ = integSimpsonLog(inf,phEmax,[e,&c,&tpf,&distCoord](double Eph)
		{
			double nPh = tpf.interpolate({{0,Eph}},&distCoord);
			double inf2 = pairThresholdPH;
			double sup2 = 2.0*Eph*e/(protonMass*cLight2);
			double integ2 = integSimpsonLog(inf2,sup2,[&](double ep)
						{
							return ep*crossSectionBetheHeitler(ep);
						},30);
			return nPh/(Eph*Eph) * integ2;
		},30);
	//return 0.5*P2(protonMass)*gsl_pow_5(cLight)/(e*e) * integ;
	return 0.5*P2(protonMass)*(cLight2*cLight2*cLight) / (e*e) * integ;
}

double pairBH(double Eee, const Particle& creator, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)
{
	double evalE = Eee*protonMass/electronMass;
	double protonDist = (evalE > creator.emin() && evalE < creator.emax()) ? 
			creator.distribution.interpolate({ { 0, evalE } }, &distCoord) : 0.0;
	
	double omega = omegaBH2(evalE,creator,tpf,distCoord,tpEmin,tpEmax);
	return 2.0*protonMass*omega*protonDist/electronMass;
}
