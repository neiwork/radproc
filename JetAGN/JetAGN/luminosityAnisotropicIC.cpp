#include "luminosityAnisotropicIC.h"

#include "lossesAnisotropicIC.h"
#include "modelParameters.h"





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

//f(x, y) en intTriple, x la de afuera, y la de adentro

double difNlum(double theta, double w, double w0, double E, double r)   //funcion a integrar
{
	// E -> energía del foton
	// w -> variable de afuera
	// w0 -> variable interna

	double b = b_theta(theta, w0, E);
	double z = E / w; 

	//defino F(z)
	double F = 1.0 + P2(z) / (2.0*(1.0 - z)) - 2.0*z / (b*(1.0 - z)) + 2.0*P2(z) / P2(b*(1.0 - z));

	double nph = blackBody(w0, r);
	double invariant = nph / w0;

	double result = (3.0*thomson / (16.0*pi)) * P2(electronMass*cLight2 / w) * invariant * F;

	return result;

}


double fLumi(double x, double theta, double y, double E,  const Particle& p, const SpaceCoord& distCoord)
{
	//double t = p.ps.current->par.T;

	double r = distCoord.dims[1]; //VER!!

	double distCreator;

	if (x < p.emin() || x > p.emax()){
		distCreator = 0.0;
	}
	else{
		distCreator = p.distribution.interpolate({ { 0, x }, { 1, r } }); //  { DIM_T, t } });// VER
	}

//	double distCreator = p.dist(u);// interpol(u, Ecreator, Ncreator, Ncreator.size() - 1);

	double result = distCreator*difNlum(theta, x, y, E, r)*(0.5*sin(theta));
	return result;
}


double luminosityAnisotropicIC(double E, Particle& particle, const SpaceCoord& distCoord)
{
	double r = distCoord.dims[1]; //VER!!

	double mass = particle.mass;

	double a = particle.emin(); //energia minima de los electrones
	double b = particle.emax();    //energia maxima de los electrones
		
	double integral =
		intTriple(E, a, b, r,
			[E](double w0, double theta){
				return cAniLum(w0, theta, E);
			},
			[E](double w0, double theta)
			{
				return dAniLum(w0, theta, E);
			},
			[E, r, &particle, &distCoord](double x, double t, double y){
				return fLumi(x, t, y, E, particle, distCoord);
			});
			
	return integral*4.0*pi*cLight*P2(E); //en erg/s/cm3

	}