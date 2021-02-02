#include "luminosityIC.h"



#include <flosses/crossSectionInel.h>
#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
//#include <fmath/interpolation.h>


#include <algorithm>


double cICemi(double u, double E, double mass, double phEmin)  //limite inferior
{
	double Erep = mass*cLight2;

	double condition = u*(u-E);

	double inf = E*P2(Erep) / (4.0*condition);

	return std::max(phEmin,inf);  //puse la condicion Ega < s*Ee/1+s
}

double dICemi(double u, double E, double phEmax)         //limite superior   
{          
	return std::min(E,phEmax);  //esta es la condicion epsilon < Ega                               
}



double fICemi2(double u, double t, double E, const Particle& creator, const SpaceCoord& distCoord,
				const ParamSpaceValues& tpf)   //funcion a integrar  u=Ee
{    
	double distCreator;
	if (u < creator.emin() || u> creator.emax()){
		distCreator = 0.0;
	}
	else{
		distCreator = creator.distribution.interpolate({ { 0, u } }, &distCoord); 
	}
	//double distCreator = creator.dist(u);// interpol(u, Ecreator, Ncreator, Ncreator.size() - 1);

	double Erep = creator.mass*cLight2;
	  

	double s = 4*u*t/P2(Erep);       //equivalente al gamma

	double r = E/(s*u*(1-E/u));   //equivalente al q

	double Nph = tpf.interpolate({ { 0, t } }, &distCoord);  

	double function = distCreator*(Nph/t)
            			*(2*r*log(r)+(1+2*r)*(1-r)+(1-r)*P2(r*s)/(2*(1+r*s)))/P2(u);
	
	double condition = s*u/(1+s);

	return function;
}


double luminosityIC(double E, const Particle& creator, const SpaceCoord& distCoord, const ParamSpaceValues& tpf,
					double phEmin, double phEmax)
{
	
	double mass = creator.mass;
	double cte  = 3.0*crossSectionThomson(creator.mass)*P2(creator.mass)*pow(cLight,5)/4;

	double integral = RungeKutta(creator.emin(), creator.emax(),
		[E,mass,phEmin](double u){
			return cICemi(u,E,mass,phEmin);
		}, 
		[E,phEmax](double u){
			return dICemi(u,E,phEmax);
		}, 
		[E,&creator,&distCoord, tpf](double u, double t){
			return fICemi2(u, t,E,creator, distCoord, tpf); 
		});    //le asigno a la variable integral el resultado de la integracion   

	double luminosity = integral*cte*P2(E);

	return luminosity;
}



double fICemi(double u, double t, double E, const Particle& creator, const SpaceCoord& distCoord, fun1 tpf)   //funcion a integrar  u=Ee
{
	double distCreator = 0.0;
	if (u >= creator.emin() && u<= creator.emax()) {
		distCreator = creator.distribution.interpolate({ { 0, u } }, &distCoord);
	}
	//double distCreator = creator.dist(u);// interpol(u, Ecreator, Ncreator, Ncreator.size() - 1);

	double Erep = creator.mass*cLight2;

	double s = 4 * u*t / P2(Erep);       //equivalente al gamma

	double r = E / (s*u*(1 - E / u));   //equivalente al q

	double Nph = tpf(t); //.interpolate({ { 0, t } }, &distCoord);

	double function = distCreator*(Nph / t)
		*(2 * r*log(r) + (1 + 2 * r)*(1 - r) + (1 - r)*P2(r*s) / (2 * (1 + r*s))) / P2(u);

	double condition = s*u / (1 + s);

	return function;
}

double luminosityIC_old(double E, const Particle& creator, const SpaceCoord& distCoord, fun1 tpf, double phEmin, double phEmax)
{

	double mass = creator.mass;
	double cte = 3.0*crossSectionThomson(creator.mass)*P2(creator.mass)*pow(cLight, 5) / 4;

	double integral = RungeKutta(creator.emin(), creator.emax(),
		[E, mass, phEmin](double u) {
		return cICemi(u, E, mass, phEmin);
	},
		[E, phEmax](double u) {
		return dICemi(u, E, phEmax);
	},
		[E, &creator, &distCoord, tpf](double u, double t) {
		return fICemi(u, t, E, creator, distCoord, tpf);
	});    //le asigno a la variable integral el resultado de la integracion   

	double luminosity = integral*cte*P2(E);

	return luminosity;
}




//**************************************************
//lo de abajo es de edu, revisar

double fICemiSy(double u, double E, const Particle& creator, const SpaceCoord& distCoord,
	double magf)   //funcion a integrar  u=Ee
{
	double distCreator;
	if (u < creator.emin() || u> creator.emax()) {
		distCreator = 0.0;
	}
	else {
		distCreator = creator.distribution.interpolate({ { 0, u } }, &distCoord);
	}
	//double distCreator = creator.dist(u);// interpol(u, Ecreator, Ncreator, Ncreator.size() - 1);

	double Erep = creator.mass*cLight2;
	double g = E / Erep;

	double Bcr = 4.4e13;
	double t = (magf / Bcr) * Erep;
	double s = 4 * u*t / P2(Erep);       //equivalente al gamma

	double r = E / (s*u*(1 - E / u));   //equivalente al q

	double NphVirtuals = magf*magf / (8 * pi) / t;

	double function = (r < 1.0 && r > 1.0 / (4 * g*g)) ? distCreator*(NphVirtuals / t)
		*(2 * r*log(r) + (1 + 2 * r)*(1 - r) + (1 - r)*P2(r*s) / (2 * (1 + r*s))) / P2(u) : 0.0;

	double condition = s*u / (1 + s);

	return function;
}
double luminositySyKN(double Eph, const Particle& creator, const SpaceCoord& distCoord, double magf)
{
	
	double mass = creator.mass;
	double cte  = 3.0*crossSectionThomson(creator.mass)*P2(creator.mass)*pow(cLight,5)/4;

	double integral = integSimpsonLog(creator.emin(), creator.emax(),
		[Eph,magf,&creator,&distCoord](double Ee){
			return fICemiSy(Ee, Eph, creator, distCoord, magf); 
		},100); 

	double luminosity = integral*cte*P2(Eph);

	return luminosity;
}

double fIC_iso(double Eg, double Ee, double Eph)
{
	double w = Eg / Ee;
	double b = 4*Eph*Ee/P2(electronRestEnergy);
	double aux = b*(1.0-w);
	return (Eg > Eph && w < b/(1+0+b)) ?
		1.0+w*w/(2.0*(1.0-w))+w/aux-2.0*w*w/P2(aux)-w*w*w/(2.0*b*P2(1.0-w))-2.0*w/aux * log(aux/w) : 0.0;
}

double luminosityIC_2(double Eg, const Particle& creator, const SpaceCoord& distCoord, const ParamSpaceValues& tpf,
					double phEmin, double phEmax)
{
	double constant = 0.75*cLight*thomson;
	double result = integSimpsonLog(creator.emin()*1.1, creator.emax()/1.1, [Eg,&creator,&distCoord,&tpf,phEmin,phEmax](double Ee)
				{
					double g = Ee / (creator.mass * cLight2);
					double ne = creator.distribution.interpolate({{0,Ee}},&distCoord);
					double w = Eg/Ee;
					//double Ephmin = max(phEmin*1.1,w/(1.0-w) * P2(electronRestEnergy)/(4*Ee));
					double Ephmin = phEmin*1.1;
					double Ephmax = min(Eg,phEmax/1.1);
					double integral = integSimpsonLog(Ephmin,Ephmax, [&distCoord,Eg,Ee,&tpf](double Eph)
					{
						double uph = tpf.interpolate({{0,Eph}},&distCoord) * Eph;
						double f = fIC_iso(Eg,Ee,Eph);
						double result1 = uph * P2(Eg/Eph) * (f > 0.0 ? f : 0.0);
						return result1;
					},50);
					return integral * ne / (g*g);
				},100);
	return constant*result;
}
