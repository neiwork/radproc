#include "lossesIC.h"

#include "crossSectionInel.h"

#include <fmath/RungeKutta.h>
#include <fparameters/parameters.h>
#include <fmath/physics.h>
//#include <gsl/gsl_sf_dilog.h>

double cIC(double u)   //limite inferior
{
	return u;
}

double dIC(double u, double E, double mass)   //limite superior     
{

	double s = 4*u*E/P2(mass*cLight2);
	return s*E/(1+s);
}

double fIC_old(double u,double t, double E, double mass, fun1 tpf)   //funcion a integrar
{
	double Erep = mass*cLight2;
	double r = t*P2(Erep) / (4 * u*E*(E - t));

	double Nterm = tpf(u);

	double result = (Nterm / u)*(t - u)*(2.0*r*log(r) +
		(1.0 + 2.0*r)*(1.0 - r) + (P2(t / (E - t))*(1.0 - r)) / (2.0*(1.0 + (t / (E - t)))));

	return  result;
}

double lossesIC_old(double E, Particle& p, fun1 tpf, double phEmin, double phEmax)

{
	//double constant  = 3*crossSectionThomson(p.mass)*P2(p.mass)*pow(cLight,5)/4;
	double constant = 3.0*thomson*p.mass*p.mass*pow(cLight, 5) / 4.0;

	double a = phEmin;     //energia minima de los fotones en erg
	double b = phEmax;     //energia maxima de los fotones en erg
	double mass = p.mass;
	double integral = RungeKutta(a, b, &cIC, [mass, E](double u) {return dIC(u, E, mass); },
		[&](double u, double t) { return fIC_old(u, t, E, mass, tpf); });    //le asigno a la variable integral el resultado de la integracion

	double de = constant*integral / P2(E);
	return de;
}



double fIC(double u,double t, double E, double mass, const ParamSpaceValues& tpf, const SpaceCoord& psc)   //funcion a integrar
{  

	double Erep = mass*cLight2;
	double r    = t*P2(Erep)/(4*u*E*(E-t));
	
	double Nterm = tpf.interpolate({ { 0, u } }, &psc	); 
	//double Nterm = tpf(u);

	double result = (Nterm/u)*(t-u)*(2.0*r*log(r)+
       	(1.0+2.0*r)*(1.0-r)+(P2(t/(E-t))*(1.0-r))/(2.0*(1.0+(t/(E-t)))));

	return  result;
}

double lossesIC(double E, Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& psc, double phEmin, double phEmax)
{
	//double constant  = 3*crossSectionThomson(p.mass)*P2(p.mass)*pow(cLight,5)/4;
	double constant = 3.0*thomson*p.mass*p.mass*pow(cLight,5)/4.0;

	double a  = phEmin;     //energia minima de los fotones en erg
	double b  = phEmax;     //energia maxima de los fotones en erg
	double mass = p.mass;
	double integral = RungeKutta(a,b,&cIC,[mass, E](double u){return dIC(u, E, mass);},
		[&](double u, double t){ return fIC(u, t, E, mass, tpf, psc); });    //le asigno a la variable integral el resultado de la integracion

	double de = constant*integral/P2(E);
	return de;
}

double lossesIC_Th(double E, Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& psc,
					double phEmin, double phEmax)
{
	double g = E / (p.mass*cLight2);
	double Uph = integSimpsonLog(phEmin, phEmax, [tpf, &psc](double Eph)
	{return tpf.interpolate({ { 0,Eph } }, &psc)*Eph; }, 100); 
	//RungeKuttaSimple(phEmin,phEmax,[tpf,&psc](double Eph)
	//{return tpf.interpolate({ {0,Eph} }, &psc)*Eph; });// , 100);
	return 4.0/3.0 * thomson * cLight * Uph * P2(electronMass/p.mass) * g*g;
}


/*
double factorQ(double e1, double e, double Ee)
{
	double gamma_e = 4.0*e*Ee/P2(electronRestEnergy);
	double q = e1/(gamma_e*(Ee-e1));
	return 2*q*log(q)+(1+2*q)*(1-q)+0.5*P2(gamma_e*q)*(1-q)/(1+gamma_e*q);
}

double dN_dtde1(double e1, double Ee, Particle& p, const ParamSpaceValues& tpf,const SpaceCoord& psc,
				double logphEmin, double phEmax)
{
	double g = Ee/electronRestEnergy;
	double Emin = 0.25*e1*electronRestEnergy/(g*(Ee-e1));
	return integSimpson(log(Emin),log(e1),[&e1,&Ee,&tpf,&psc](double loge)
			{
				double e = exp(loge);
				return tpf.interpolate({{0,e}},&psc) * factorQ(e1,e,Ee);
			},40);
}



double fKN(double b)
{
	double g = (0.5*b + 6.0 + 6.0/b) * log(1.0+b);
	g = g - (11.0/12.0 * b*b*b + 6.0*b*b + 9.0*b + 4.0) / P2(1.0+b);
	g = g - 2.0 + 2.0 * gsl_sf_dilog(-b);
	return (b > 1e-3 ? 9.0*g/(b*b*b) : 1.0);
}

double lossesIC(double E, Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& psc,
					double phEmin, double phEmax)
{
	double g = E / (p.mass*cLight2);
	double Uph = integSimpsonLog(phEmin,phEmax,[tpf,&psc](double Eph)
					{return tpf.interpolate({{0,Eph}},&psc)*Eph;},100);
	double lossesTh = 4.0/3.0 * thomson * cLight * Uph * P2(electronMass/p.mass) * g*g;
	return (lossesTh/Uph) * integSimpsonLog(phEmin,phEmax, [&tpf,&psc,g,&p] (double Eph)
				{
					double nPh = tpf.interpolate({{0,Eph}},&psc);
					double b = 4.0 * g * (Eph/(p.mass*cLight2));
					return fKN(b)*nPh*Eph;
				},100);
}
*/