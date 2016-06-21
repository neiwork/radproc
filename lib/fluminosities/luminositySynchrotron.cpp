#include "luminositySynchrotron.h"


#include "opticalDepthSSA.h"
#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\physics.h>



double fSyn(double x, double E, const Particle& creator, const ParamSpaceValues& magf, const SpaceCoord& psc)         //funcion a integrar   x=Ee; L=L(Ega)
{
	//double r = creator.ps.current->par.R;
	//double t = creator.ps.current->par.T;

	const double magneticField = magf.get(psc);
	double distCreator;
	if (x < creator.emin() || x> creator.emax()){
		distCreator = 0.0;
	}
	else{
		distCreator = creator.distribution.interpolate({ { 0, x } }, &psc	); 
	}
//	double distCreator = creator.dist(x);// interpol(x, Ecreator, Ncreator, Ncreator.size() - 1);

	double cte = pow(3.0, 0.5)*P3(electronCharge)*magneticField / (planck*creator.mass*cLight2);

	double Echar = 3 * electronCharge*planck*magneticField*P2(x) / (4 * pi*P3(creator.mass)*cLight*P2(cLight2));
	
	double aux = E/Echar;  //aca el aux es el x real

	double result = cte*1.85*distCreator*pow(aux,(1.0/3.0))*exp(-aux);  

	return result;
}


double luminositySynchrotron(double E, const Particle& c, const SpaceCoord& psc, const ParamSpaceValues& magf)
{
	
	double integralS = RungeKuttaSimple(c.emin(), c.emax(), [&](double e){
		return fSyn(e, E, c, magf, psc);
	});

	double luminosityS = integralS*E; //multiplico por E asi obtengo luminosidad

	if (luminosityS > 0.0){ return luminosityS; }
	else { return 0.0; }

}




//double luminositySynchrotron_conSSA(double E, const Particle& creator)
//{
//
//	double Emax = creator.emax();
//	double Emin = creator.emin();
//
//	double tau = opticalDepthSSA(E, creator.mass, Emin, Emax, creator);  //E=Eph
//
//	double factorSSA = (1.0-exp(-tau))/tau;
//
//	if (factorSSA > 1.0 || factorSSA == 0.0) //1.0e-3)  //lo cambio solo en el caso que interesa
//	{ factorSSA = 1.0; }	
//	
//	//double integralS = RungeKuttaSimple(Emin, Emax, bind(fSyn,_1,E,creator));
//
//	double integralS = RungeKuttaSimple(creator.emin(), creator.emax(),   //RungeKuttaSimple(double a, double b, fun1 f)
//		[E, &creator](double x){
//		return fSyn(x, E, creator);  //double fSyn(double x, double E, const Particle& creator)
//	});
//
//
//	double luminosityS = factorSSA*integralS*E; //multiplico por E asi obtengo luminosidad
//
//	if (luminosityS > 0.0){return luminosityS ; }
//	else {return 0.0;} 
//	 
//}


////////////////////////////////////////////////////////////////////////////////////////////////////
/*
double fSynS(double x, double E, const Particle& creator)         //funcion a integrar   x=Ee; L=L(Ega)
{
	double distCreator = creator.dist(x);//interpol(x,Ecreator,Ncreator,Ncreator.size()-1);

	double cte	= pow(3.0,0.5)*P3(electronCharge)*magneticField/(planck*creator.mass*cLight2);

	double Echar = 3*electronCharge*planck*magneticField*P2(x)/(4*pi*P3(creator.mass)*cLight*P2(cLight2));
	
	double tau = E/x;
	double aux = E/(Echar*(1-tau));  //aca el aux es el x real

	double result = cte*1.85*(1-tau)*distCreator*pow(aux,(1.0/3.0))*exp(-aux);  

	return result;    // esta condicion la puse en el limite inferior tau<1 ? result : 0.0;
}



double luminositySynchrotronSec(double E, const Particle& c)
{
	double Emax = c.emax();
	double Emin = c.emin();
	double inf  = std::max(Emin,E);   //esto lo agrego asi le saco la condicion sobre tau < 1

	double tau = opticalDepthSSA(E, c.mass, inf, Emax, c);  //E=Eph

	double factorSSA = 1.0;

	if (tau > 1.0e-3)  //lo cambio solo en el caso que interesa
	{ factorSSA = (1.0-exp(-tau))/tau;}

	double integralS = RungeKuttaSimple(Emin, Emax, [&E, &c](double e){
		return fSynS(e, E, c);
	});

	double luminosityS = factorSSA*integralS*E*volume; //multiplico por E asi obtengo luminosidad
	                                 //divido por E asi obtengo emisividad y no luminosidad

	return luminosityS ; 

}

*/