#include "muonInj.h"


#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
#include <fmath/interpolation.h>
//#include <algorithm>


double fQ_L(double Epi, double E, const Particle& p, const Particle& c, const SpaceCoord& psc)  //funcion a integrar variable Epi
{
	
	double distPion = c.distribution.interpolate({ { 0, Epi } }, &psc);  //c.dist(Epi);

	double r = P2(p.mass/chargedPionMass);

	double x = E/Epi;

	double decayTime = chargedPionMeanLife*Epi/(chargedPionMass*cLight2);

	double Q_L;

//	if(x > r && x < 1)	{                //esta condición es por la función de Heaviside
										 //la condición de x es para evitar errores numericos
		Q_L = distPion*r*(1-x)/(Epi*x*P2(1-r)*decayTime);   
//	}
//	else { Q_L = 0.0;	}
	
	return Q_L;    //esto es muL-  y muR+
}

double fQ_R(double Epi, double E, const Particle& p, const Particle& c, const SpaceCoord& psc)         //funcion a integrar variable Epi
{

	double distPion = c.distribution.interpolate({ { 0, Epi } }, &psc); //dist(Epi);

	double r = P2(p.mass/chargedPionMass);

	double x = E/Epi;

	double decayTime = chargedPionMeanLife*Epi/(chargedPionMass*cLight2);

	double Q_R;  

//	if(x > r && x < 1)	{                //esta condición es por la función de Heaviside
										  //la condición de x es para evitar errores numericos
		Q_R = distPion*(x-r)/(Epi*x*P2(1-r)*decayTime);
//	}
//	else { Q_R = 0.0;	}
	
	return Q_R;   //esto es muR- y muL+
}


double muonInj(double E, const Particle& p, const Particle& c, const SpaceCoord& psc)  //el proton es el creator del pion
{        
	
	std::string pName = p.id;
	double rpi = P2(muonMass/chargedPionMass);
	double injection=0.0;
	double sup = std::min(c.emax(),E/rpi);  // transformo la condicion de la heaviside en un limite superior

	fun1 rk_fQl = [&p,&c,E,&psc](double e){
		return fQ_L(e, E, p, c, psc);
	};

	fun1 rk_fQr = [&p, &c, E,&psc](double e){
		return fQ_R(e, E, p, c,psc);
	};

	if (p.id == "ntMuon") {
		injection = RungeKuttaSimple(E, sup, rk_fQl);
		injection += RungeKuttaSimple(E, sup, rk_fQr);//*2; 
		//el dos es por las dos particulas para cada caso
	}
	else if (p.id == "muon_L_minus" && pName == "muon_R_plus") {  //muon_L- y muon_R+
		injection = RungeKuttaSimple(E, sup, rk_fQl);
	}
	else if (p.id == "muon_R_minus" && pName == "muon_L_plus") {  //muon_L+ y muon_R-
		injection = RungeKuttaSimple(E, sup, rk_fQr);
	}
	return injection;

}

	/*switch (particleName)	{
		case PT_muon_L_minus:  //muon_L- y muon_R+
		case PT_muon_R_plus:	
			injection = RungeKuttaSimple(E, sup, rk_fQl);
			break;
		case PT_muon_R_minus:  //muon_L+ y muon_R-
		case PT_muon_L_plus:
			injection = RungeKuttaSimple(E, sup, rk_fQr);
			break;
		case PT_muon:   
			injection = RungeKuttaSimple(E, sup, rk_fQl);
			injection += RungeKuttaSimple(E, sup, rk_fQr);//*2; 
			break;    //el dos es por las dos particulas para cada caso
	}*/


double fQ_L2(double Epi, double Emu, const Particle& p)  //funcion a integrar variable Epi
{
	double r = P2(p.mass/chargedPionMass);
	double x = Emu/Epi;
	return (x > r) ? r*(1.0-x)/(Epi*x*P2(1.0-r)) : 0.0;
}

double fQ_R2(double Epi, double Emu, const Particle& p)
{
	double r = P2(p.mass/chargedPionMass);
	double x = Emu/Epi;
	return (x > r) ? (x-r)/(Epi*x*P2(1.0-r)) : 0.0;
}

double muonInjNew(double Emu, const Particle& p, const Particle& c, const SpaceCoord& psc)  //el proton es el creator del pion
{        
	/*double injection = integSimpson(log(Emu),log(c.emax()),[Emu,&p,&c,&psc](double logEpi)
						{
							double Epi = exp(logEpi);
							double tDecay = chargedPionMeanLife*(Epi/(chargedPionMass*cLight2));
							double nPi = c.distribution.interpolate({{0,Epi}},&psc);
							return Epi*nPi/tDecay * (fQ_R2(Epi,Emu,p)+fQ_L2(Epi,Emu,p));
						},50);*/
	double injection = integSimpsonLog(Emu,c.emax(),[Emu,&p,&c,&psc](double Epi)
						{
							double tDecay = chargedPionMeanLife*(Epi/(chargedPionMass*cLight2));
							double nPi = c.distribution.interpolate({{0,Epi}},&psc);
							return nPi/tDecay * (fQ_R2(Epi,Emu,p)+fQ_L2(Epi,Emu,p));
						},50);
	return injection;
}