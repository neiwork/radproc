#include "luminosityAnisotropicIC.h"

//#include "lossesAnisotropicIC.h"
#include "modelParameters.h"

#include "targetFields.h"

#include <fparameters\SpaceCoord.h>
#include <fparameters\Dimension.h>
#include <algorithm>

#include "modelParameters.h"
#include <fmath\physics.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>


//u=epsilon
//b = beta
//eGamma electron Lorentz factor, E/mc^2
//E = Egamma = es

double arg(double es, double u, double g, double eBeta, double k,
	double eta, double eta_e, double eta_obs, double phi)
{
	double num = u*g*(1.0 - eBeta*eta*eta_e - eBeta*sqrt((1.0 - P2(eta))*(1.0 - P2(eta_e)))*cos(phi));
	double den = (u*(1.0 - k) + g*(1.0 - eBeta*eta_obs*eta_e));
	double res = es - num / den;
	return res;
}

double r_eta(double eta, double z, double bG)
{
	double r = -sqrt(P2(eta*bG) - P2(eta) - P2(bG) + 1)*z / (eta + bG);
	return (r > 0) ? r : -r;

	//double r = z*sqrt(P2(eta*bG) +2.0*eta*bG)/ (eta + bG);
	//return r; // (r > 0) ? r : -r;

}



double eta_ph_z(double r, double z, double bG)
{
	double x = sqrt(P2(r) + P2(z));
	double eta_ph = (z - bG*x) / (x - bG*z);
	return eta_ph;

}
double mu_i(double eta_e, double eta, double phi)
{
	double mu = eta*eta_e + sqrt((1.0 - P2(eta))*(1.0 - P2(eta_e)))*cos(phi);
	return mu;
}

double f_sum(double eta_e, double eta, double eta_obs, double u,
	double phi, double g, double eBeta, double k, double es)
{
	double mu = mu_i(eta_e, eta, phi);
	double mu_obs = eta_obs*eta_e;

	double dude = eta - eta_e*cos(phi)*sqrt((1 - P2(eta)) / (1 - P2(eta_e)));

	double bnn = eBeta*eta_obs*eta_e;

	double num = u*g*eBeta*dude*(u*(1.0 - k) + g*(1.0 - bnn))
		- u*P2(g)*eBeta*eta_obs*(1.0 - eBeta*mu);

	double den = P2(u*(1.0 - k) + g*(1.0 - bnn));

	double dfdn = abs(num / den);

	double factor = (1.0 - eBeta*mu)*dfdn;  

	double f1 = u*(1.0 - eBeta*mu) / (es*(1.0 - eBeta*mu_obs));

	double kp = 1.0 - (1.0 - k) / (P2(g)*(1.0 - eBeta*mu_obs)*(1.0 - eBeta*mu));
	double sin_xi = sqrt(1.0 - P2(kp));

	return (f1 + 1.0/f1 - 1.0+P2(kp)) / factor;
}



double luminosityAnisotropicIC(double E, Particle& particle, double z,
	double gamma, fun3 tpf, double starT, const SpaceCoord& psc)
{
	static const double R_d = GlobalConfig.get<double>("R_d") *pc;
	static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree

	double mass = particle.mass;
	double Erest = (mass*cLight2);

	double es = E/Erest; //es la energia del foton gamma
	
	double bG = beta(gamma);
	double eta_obs = (cos(inc) - bG) / (1.0 - bG*cos(inc)); //cos(inc) = eta_obs*
	
	double theta = boltzmann*starT / Erest;

	double g_min = particle.emin()/Erest;
	double g_max = particle.emax() / Erest;
	int n_g = 30;
	double g_int = pow((g_max / g_min), (1.0 / n_g));
	
	double phi_min = 0.0;
	double phi_max = 2.0*pi;
	int nPhi = 10;
	double dphi = (phi_max - phi_min) / nPhi;  //este es phi_ph

	double eta_min = eta_ph_z(R_d, z, bG); //eta_min con Rmax, y viceversa; deta/dR < 0
	double eta_max = eta_ph_z(1.0e11, z, bG); //1.0e11 -> generalizar
	int nEta = 10;
	double deta = (eta_max - eta_min) / nEta;

	double suma = 0.0;

	double g = g_min;
	// [hook]
	// { the only sincronization point is the accumulation in suma }
	// { we use the reduction directive to specify this: }
#pragma omp parallel for reduction(+:suma)

	for (int i_g = 0; i_g < n_g; ++i_g)
	{
		double eBeta = beta(g);  //g = factor lorentz electrones
		double dg = g*(g_int - 1.0);

		double Ne = particle.distribution.interpolate({ { 0, g*Erest } }, &psc) / (4.0*pi);

		double f_g = Ne*Erest / P2(g); //paso N(E) -> N(g) 

		double phi = phi_min;
		for (int i_phi = 0; i_phi < nPhi+1; ++i_phi)
		{

			double eta = eta_min;  //este eta es en realidad eta_ph
			for (int i_eta = 0; i_eta < nEta+1; ++i_eta)
			{

				double u = 2.7*theta / (gamma*(1.0 + bG*eta));  //<eps>
				
				double k = eta*eta_obs+sqrt((1.0 - P2(eta))*(1.0 - P2(eta_obs)))*
							cos(phi-pi/2.0);
				                   
				double cc = 1.0 - P2(eta);

				double alpha = P2(u*g*eBeta*cos(phi))*cc;
				double rho = es*u*(1.0 - k) + g*(es - u);
				double pseda = g*eBeta*(u*eta - es*eta_obs);
				
				double a =  alpha + P2(pseda);
				double b = 2.0*rho*pseda;
				double c = P2(rho) - alpha; 
			
				double q;
				double basc = P2(b) - 4.0*a*c;
				double basc2 = P2(alpha) + alpha*(P2(pseda)-P2(rho));
				
				
				if (basc > 0.0){

					q = -0.5*(b + sgn(b)*sqrt((basc)));

					double eta_e1 = q / a;
					double eta_e2 = c / q;

					//double eta_e1 = (-pseda*rho + sqrt(abs(basc2))) / (2.0*(P2(pseda+alpha)));
					//double eta_e2 = (-pseda*rho - sqrt(abs(basc2))) / (2.0*(P2(pseda + alpha)));

					double f1 = 0.0;
					double f2 = 0.0;

					if (abs(eta_e1) <= 1.0)
					{
						f1 = f_sum(eta_e1, eta, eta_obs, u, phi, g, eBeta, k, es);
						double delta = arg(es, u, g, eBeta, k, eta, eta_e1, eta_obs, phi);
					}
					if( abs(eta_e2) <= 1.0){
						f2 = f_sum(eta_e2, eta, eta_obs, u, phi, g, eBeta, k, es);
						double delta = arg(es, u, g, eBeta, k, eta, eta_e2, eta_obs, phi);
					}
					
					double f_t = P3(theta / (1.0 + bG*eta));

					double integral = f_g*f_t*P2(es / u)*(f1 + f2);

					double r = r_eta(eta, z, bG);
					double r2 = r_eta(eta-deta, z, bG);
					double dr = abs(r-r2);
					//double x = sqrt(P2(r) + P2(z));
					//double dr = deta*x*P2(gamma*(x-bG*z))/(z*r);
					
					double nph_norm = dr*tpf(u*Erest, z, r); // (u*Erest); 
															// [nph_norm]= cm^-3, en el sist del disco

					suma += integral*nph_norm*dg*dphi*deta;
				}
				eta = eta + deta; // *eta_int;
			}  //cierro eta

			phi = phi + dphi;
		} //cierro phi

		g = g*g_int;
	} //cierro r

	double cte = 3.0*thomson*cLight*2.404/(8.0*P3(gamma));

	return suma*cte;// *P2(E); // *P2(E) para que quede en erg/s
							
							  
}


#include "luminosityAnisotropicIC.h"

#include "write.h"
//#include "lossesAnisotropicIC.h"
#include "modelParameters.h"

#include "targetFields.h"

#include <fparameters\SpaceCoord.h>
#include <fparameters\Dimension.h>
#include <algorithm>

#include "modelParameters.h"
#include <fmath\physics.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>


//u=epsilon
//b = beta
//eGamma electron Lorentz factor, E/mc^2
//E = Egamma = es

double arg(double es, double u, double g, double eBeta, double k,
	double eta, double eta_e, double eta_obs, double phi)
{
	double num = u*g*(1.0 - eBeta*eta*eta_e - eBeta*sqrt((1.0 - P2(eta))*(1.0 - P2(eta_e)))*cos(phi));
	double den = (u*(1.0 - k) + g*(1.0 - eBeta*eta_obs*eta_e));
	double res = es - num / den;
	return res;
}

double r_eta(double eta, double z, double bG)
{
	double r = -sqrt(P2(eta*bG) - P2(eta) - P2(bG) + 1)*z / (eta + bG);
	return (r > 0) ? r : -r;

	//double r = z*sqrt(P2(eta*bG) +2.0*eta*bG)/ (eta + bG);
	//return r; // (r > 0) ? r : -r;

}



double eta_ph_z(double r, double z, double bG)
{
	double x = sqrt(P2(r) + P2(z));
	double eta_ph = (z - bG*x) / (x - bG*z);
	return eta_ph;

}
double mu_i(double eta_e, double eta, double phi)
{
	double mu = eta*eta_e + sqrt(abs(1.0 - P2(eta))*abs(1.0 - P2(eta_e)))*cos(phi);
	return mu;
}

double f_sum(double eta_e, double eta, double eta_obs, double u,
	double phi, double g, double eBeta, double k, double es)
{
	double mu = mu_i(eta_e, eta, phi);
	double mu_obs = eta_obs*eta_e;

	double dude = eta - eta_e*cos(phi)*sqrt(abs(1 - P2(eta)) / abs(1 - P2(eta_e)));

	double bnn = eBeta*eta_obs*eta_e;

	double num = u*g*eBeta*dude*(u*(1.0 - k) + g*(1.0 - bnn))
		- u*P2(g)*eBeta*eta_obs*(1.0 - eBeta*mu);

	double den = P2(u*(1.0 - k) + g*(1.0 - bnn));

	double dfdn = abs(num / den);

	double factor = (1.0 - eBeta*mu)*dfdn;  

	double f1 = u*(1.0 - eBeta*mu) / (es*(1.0 - eBeta*mu_obs));

	double kp = 1.0 - (1.0 - k) / (P2(g)*(1.0 - eBeta*mu_obs)*(1.0 - eBeta*mu));
	double sin_xi = sqrt(1.0 - P2(kp));

	return (f1 + 1.0/f1 - 1.0+P2(kp)) / factor;
}



double luminosityAnisotropicIC(double E, Particle& particle, double z,
	double gamma, fun3 tpf, double starT, const SpaceCoord& psc)
{
	static const double R_d = GlobalConfig.get<double>("R_d") *pc;
	static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree

	double mass = particle.mass;
	double Erest = (mass*cLight2);

	double es = E / Erest; //es la energia del foton gamma

	double bG = beta(gamma);
	double eta_obs = (cos(inc) - bG) / (1.0 - bG*cos(inc)); //cos(inc) = eta_obs*

	double theta = boltzmann*starT / Erest;

	double g_min = particle.emin() / Erest;
	double g_max = particle.emax() / Erest;
	int n_g = 30;
	double g_int = pow((g_max / g_min), (1.0 / n_g));

	double phi_min = 0.0;
	double phi_max = 2.0*pi;
	int nPhi = 10;
	double dphi = (phi_max - phi_min) / nPhi;  //este es phi_ph

	double eta_min = eta_ph_z(R_d, z, bG); //eta_min con Rmax, y viceversa; deta/dR < 0
	double eta_max = eta_ph_z(pc, z, bG); //1.0e11 -> generalizar
	int nEta = 10;
	double deta = (eta_max - eta_min) / nEta;

	double suma = 0.0;

	double g = g_min;
	// [hook]
	// { the only sincronization point is the accumulation in suma }
	// { we use the reduction directive to specify this: }
#pragma omp parallel for reduction(+:suma)

	for (int i_g = 0; i_g < n_g; ++i_g)
	{
		double eBeta = beta(g);  //g = factor lorentz electrones
		double dg = g*(g_int - 1.0);

		double Ne = particle.distribution.interpolate({ { 0, g*Erest } }, &psc) / (4.0*pi);

		double f_g = Ne*Erest / P2(g); //paso N(E) -> N(g) 

		double phi = phi_min;
		for (int i_phi = 0; i_phi < nPhi + 1; ++i_phi)
		{
			double eta = eta_min;  //este eta es en realidad eta_ph
			for (int i_eta = 0; i_eta < nEta + 1; ++i_eta)
			{

				double u = 2.7*theta / (gamma*(1.0 + bG*eta));  //<eps>

				double k = eta*eta_obs + sqrt(abs(1.0 - P2(eta))*abs(1.0 - P2(eta_obs)))*
					cos(phi - pi / 2.0);
				
				double c = P2(es*u*(1.0-k)+g*(es-u))-P2(g*eBeta*cos(phi))*abs(1.0 - P2(eta));
				double b = -2.0*g*eBeta*(es*eta_obs - u*eta)*(es*u*(1.0 - k) + g*(es - u));
				double a = P2(g*eBeta*(es*eta_obs-u*eta))+ P2(g*eBeta*cos(phi))*abs(1.0-P2(eta));
								
				double basc1 = P2(es*u*(1.0 - k) + g*(es - u));
				double basc2 = abs(basc1/a);
				
				double eta_e1 = 0.0;
				double eta_e2 = 0.0;
				double f1 = 0.0;
				double f2 = 0.0;

				if (abs(basc2 - 1.0) < 1.0e-5) {
					eta_e1 = -b / (2.0*a); //basc = 0
					eta_e2 = 0.0;					
					f1 = f_sum(eta_e1, eta, eta_obs, u, phi, g, eBeta, k, es);
					
				}
				else if (a - basc1 > 0.0) {
					eta_e1 = (-b+sqrt(P2(b)-4.0*a*c))/(2.0*a);
					eta_e2 = (-b - sqrt(P2(b) - 4.0*a*c)) / (2.0*a);

					f1 = f_sum(eta_e1, eta, eta_obs, u, phi, g, eBeta, k, es);
					f2 = f_sum(eta_e2, eta, eta_obs, u, phi, g, eBeta, k, es);
				}

				if (abs(eta_e1 - 1.0) < 1.0e-3 || abs(eta_e1) >= 1.0 )
				{
					eta_e1 = 1.0;
					f1 = 0.0;
				}

				if (abs(eta_e2 - 1.0) < 1.0e-3 || abs(eta_e2) >= 1.0)
				{
					eta_e2 = 1.0;
					f2 = 0.0;
				}

				double f_t = P3(theta / (1.0 + bG*eta));

				double integral = f_g*f_t*P2(es / u)*(f1 + f2);

				double r = r_eta(eta, z, bG);
				double r2 = r_eta(eta - deta, z, bG);
				double dr = abs(r - r2);

				double nph_norm = dr*tpf(u*Erest, z, r); // (u*Erest); 
														// [nph_norm]= cm^-3, en el sist del disco
				if (integral > 0.0) {
					suma += integral*nph_norm*dg*dphi*deta;
				}
				else {

					std::cout << eta_e1 << '\t' << eta_e2 << std::endl;
					std::cout << f1 << '\t' << f2 << std::endl;
				}
				eta = eta + deta; // *eta_int;
			}  //cierro eta

			phi = phi + dphi;
		} //cierro phi

		g = g*g_int;
	} //cierro r

	double cte = 3.0*thomson*cLight*2.404 / (8.0*P3(gamma));

	return suma*cte;// *P2(E); // *P2(E) para que quede en erg/s


}





/*double luminosityAnisotropicIC(double E, Particle& particle, double z,
	double gamma, fun3 tpf, double starT, const SpaceCoord& psc)
{
	static const double R_d = GlobalConfig.get<double>("R_d") *pc;
	static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree

	double mass = particle.mass;
	double Erest = (mass*cLight2);

	double es = E/Erest; //es la energia del foton gamma
	
	double bG = beta(gamma);
	double eta_obs = (cos(inc) - bG) / (1.0 - bG*cos(inc)); //cos(inc) = eta_obs*
	
	double theta = boltzmann*starT / Erest;

	double cond1 = es / (gamma*(1.0 - bG*(cos(inc))));
	double g_min = std::max(particle.emin()/Erest,cond1);
	double g_max = particle.emax() / Erest;
	int n_g = 50;
	double g_int = pow((g_max / g_min), (1.0 / n_g));
	
	double phi_min = 0.0;
	double phi_max = 2.0*pi;
	int nPhi = 20;
	double dphi = (phi_max - phi_min) / nPhi;  //este es phi_ph

	double eta_ph_min = eta_ph_z(R_d, z, bG); //eta_min con Rmax, y viceversa; deta/dR < 0
	double eta_ph_max = eta_ph_z(1.0e11, z, bG); //1.0e11 -> generalizar
	int nEta = 20;
	double detaph = (eta_ph_max - eta_ph_min) / nEta;

	double eta_e_min = -1.0; 
	double eta_e_max = 1.0; 
	int nEta_e = 50;
	double detae = (eta_e_max - eta_e_min) / nEta_e;



	double suma = 0.0;

	double g = g_min;
	// [hook]
	// { the only sincronization point is the accumulation in suma }
	// { we use the reduction directive to specify this: }
#pragma omp parallel for reduction(+:suma)

	for (int i_g = 0; i_g < n_g; ++i_g)
	{
		double eBeta = beta(g);  //g = factor lorentz electrones
		double dg = g*(g_int - 1.0);

		double Ne = particle.distribution.interpolate({ { 0, g*Erest } }, &psc) / (4.0*pi);

		double f_g = Ne*Erest / P2(g); //paso N(E) -> N(g) 

		double phi = phi_min;
		for (int i_phi = 0; i_phi < nPhi+1; ++i_phi)
		{

			double eta_ph = eta_ph_min;  //este eta es en realidad eta_ph
			for (int i_eta_ph = 0; i_eta_ph < nEta+1; ++i_eta_ph)
			{

				double u = 2.7*theta / (gamma*(1.0 + bG*eta_ph));  //<eps>
				
				double cond2 = u*es*(1.0 - eta_obs*eta_ph);

				double k = eta_ph*eta_obs+sqrt((1.0 - P2(eta_ph))*(1.0 - P2(eta_obs)))*
							cos(phi-pi/2.0);

				double r = r_eta(eta_ph, z, bG);
				double r2 = r_eta(eta_ph-detaph, z, bG);
				double dr = abs(r-r2);
					
				double nph_norm = dr*tpf(u*Erest, z, r); // (u*Erest); 
														// [nph_norm]= cm^-3, en el sist del disco

				double f_ph = nph_norm*P3(theta / (1.0 + bG*eta_ph));
					
				double eta_e = eta_e_min;  
				for (int i_eta_e = 0; i_eta_e < nEta_e + 1; ++i_eta_e)
				{
					
					double delta = arg(es, u, g, eBeta, k, eta_ph, eta_e, eta_obs, phi);
					//std::cout << delta << '\t' << eta_e << '\t' << log10(E/1.6e-12) << std::endl;

					double mu_obs = eta_obs*eta_e;

					double mu = mu_i(eta_e, eta_ph, phi);

					double prueba = es*u*(1.0 - mu*mu_obs);

					if (abs(delta) < 1.0e-15) {

						

						double f1 = u*(1.0 - eBeta*mu) / (es*(1.0 - eBeta*mu_obs));

						double kp = 1.0 - (1.0 - k) / (P2(g)*(1.0 - eBeta*mu_obs)*(1.0 - eBeta*mu));
						double sin_xi = sqrt(1.0 - P2(kp));

						return (f1 + 1.0 / f1 - 1.0 + P2(kp));

						double integral = f_g*(f_ph / (1.0 - eBeta*mu))*P2(es / u);


						suma += integral*nph_norm*dg*dphi*detaph*detae;
					}

					eta_e = eta_e + detae;
				}
				eta_ph = eta_ph + detaph; // *eta_int;
			}  //cierro eta

			phi = phi + dphi;
		} //cierro phi

		g = g*g_int;
	} //cierro r

	double cte = 3.0*thomson*cLight*2.404/(8.0*P3(gamma));

	return suma*cte;// *P2(E); // *P2(E) para que quede en erg/s*/