#include "luminosityAnisotropicIC.h"

#include "lossesAnisotropicIC.h"





double cAniLum(double w0, double theta, double E)   //limite inferior
{
	double s = b_theta(theta, w0, E);  //reemplazo el Gamma = 4*u*E / P2(mass*cLight2), para incuir la dep theta
	double inf =  s*E / (1 + s);
	
	inf =std:: max(inf, targetPhotonEmin);
	
	return inf;   ////puse la condicion Ega < s*Ee/1+s, recordar que ahora E = Ega
}

double dAniLum(double w0, double theta, double E)   //limite superior     
{
	return E;  //esta es la condicion epsilon < Ega	= E
}


double fLumi(double x, double theta, double y, double E, double r, const Particle& p)
{
	//p.ps.dimensions[2]->values[])
	double distCreator = p.distribution.interpolate({ x, r, p.ps[2][0] });//p.ps[2].par.T;// VER

//	double distCreator = p.dist(u);// interpol(u, Ecreator, Ncreator, Ncreator.size() - 1);

	return distCreator*difN(theta, x, y, E, r)*(0.5*sin(theta));
}


double luminosityAnisotropicIC(double E, Particle& particle, double r)
{
	double mass = particle.mass;

	double a = pow(10.0,particle.logEmin)*1.6e-12;      //energia minima de los electrones
	double b = pow(10.0, particle.logEmax)*1.6e-12;    //energia maxima de los electrones
		
	double integral =
		intTriple(E, a, b, r,
		[E](double w0, double theta){
		return cAniLum(w0, theta, E); },
			[E](double w0, double theta){
			return dAniLum(w0, theta, E); },
				[E, r, &particle](double x, double t, double y){
				return fLumi(x, t, y, E, r, particle); });
			
	return integral*4.0*pi*cLight*P2(E); //en erg/s

	}



