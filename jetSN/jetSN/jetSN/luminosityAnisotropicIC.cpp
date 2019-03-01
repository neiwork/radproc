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

double KNcrossSection(double eval)
{
	double cs = 3.0*thomson*(log(2.0*eval) + 0.5)/(8.0*eval);
	return cs;
}

double eta_ph_z(double r, double z, double bG)
{
	double x = sqrt(P2(r) + P2(z));
	double eta_ph = (z - bG*x) / (x - bG*z);
	return eta_ph;

}

double r_eta(double eta, double z, double bG)
{
	double r = -sqrt(P2(eta*bG) - P2(eta) - P2(bG) + 1)*z / (eta + bG);
	return (r > 0) ? r : -r;

	//double r = z*sqrt(P2(eta*bG) +2.0*eta*bG)/ (eta + bG);
	//return r; // (r > 0) ? r : -r;

}

double ne_eval(Particle& particle, double eval_e, const SpaceCoord& psc)
{
	double Erest = particle.mass*cLight2;
	double Ne = 0.0;

	double emin = particle.emin();
	double emax = particle.emax();

	if (eval_e <= emax && eval_e >= emin) {
		Ne = Erest*particle.distribution.interpolate({ { 0, eval_e } }, &psc) / (2.0*pi);
	}

	return Ne;
}


double luminosityAnisotropicIC(double E, Particle& particle, double z,
	double gamma, fun3 tpf, double starT, const SpaceCoord& psc)
{
	static const double R_d = GlobalConfig.get<double>("R_d") *pc;
	static const double h_d = GlobalConfig.get<double>("h_d") *pc;
	static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree
	static const double logEmax = GlobalConfig.get<double>
		("model.particle.photon.dim.energy.max", 1.0);
	
	double photonEmax = pow(10.0,logEmax)*1.6e-12;

	double D = computeDlorentz(gamma);

	double Erest = (particle.mass*cLight2);

	const double es = E / Erest; //es la energia del foton gamma

	const double bG = beta(gamma);
	const double eta_obs = (cos(inc) - bG) / (1.0 - bG*cos(inc)); //cos(inc) = eta_obs*

	const double theta = boltzmann*starT / Erest;

	double eta_min = eta_ph_z(R_d, z, bG); //eta_min con Rmax, y viceversa; deta/dR < 0
	double eta_max = eta_ph_z(1.0e15, z, bG); //1.0e11 -> generalizar
	int nEta = 20;
	double deta = (eta_max - eta_min) / nEta;

	double phi_min = 0.0;
	double phi_max = 2.0*pi;
	int nPhi = 20;
	double dphi = (phi_max - phi_min) / nPhi;

	double suma = 0.0;

	double eta = eta_min;
	// [hook]
	// { the only sincronization point is the accumulation in suma }
	// { we use the reduction directive to specify this: }
#pragma omp parallel for reduction(+:suma)

	for (int i_eta = 0; i_eta < nEta + 1; ++i_eta)
	{
		double phi = phi_min;
		for (int i_phi = 0; i_phi < nPhi + 1; ++i_phi)
		{
			double mu = eta*eta_obs + sqrt(abs(1.0 - P2(eta))*abs(1.0 - P2(eta_obs)))*cos(phi);
		//    integral 1, regimen thomson

			double g_min = particle.emin() / Erest;
			g_min = std::max(g_min, es);
			double g_max = particle.emax() / Erest;
			int n_g = 80;
			double g_int = pow((g_max / g_min), (1.0 / n_g));

			double suma_1 = 0.0;
			double g = g_min;
			for (int i_g = 0; i_g < n_g; ++i_g)
			{
				double eBeta = beta(g);  //g = factor lorentz electrones
				double dg = g*(g_int - 1.0);

				double eval_e = g*Erest;
				double Ne = ne_eval(particle, eval_e, psc);

				double f_g = Ne / P2(g); 
	
				//double eval_ph = gamma*es*(1.0 + bG*eta) / (P2(g)*(1.0 - eta*eta_obs));
				double eval_ph = gamma*es*(1.0 + bG*eta) / (P2(g)*(1.0 - eBeta*mu));

				double nph = tpf(eval_ph, gamma, eta);// nph_ICani2(eval_ph, gamma, eta, 'IR') ; //falta dr

				suma_1 += cLight*thomson*f_g*nph*dg;


				g = g*g_int;
			}

		//    integral 2, sobre eps, regimen KN
			double suma_2 = 0.0;

			//double e_min = std::max(1.0, 1.0 / (es*(1.0 - eta*eta_obs)));
			double e_min = std::max(1.0, 1.0 / (es*(1.0 - mu)));
			
			double e_max = photonEmax / Erest;
			int n_e = 80;
			double e_int = pow((e_max / e_min), (1.0 / n_e));

			double eps = e_min;
			for (int i_e = 0; i_e < n_e; ++i_e)
			{

				double de = eps*(e_int - 1.0);

				double eval_ph = gamma*eps*(1.0 + bG*eta);
				double nph = tpf(eval_ph, gamma, eta);
				// nph_ICani2(eval_ph, gamma, eta, 'IR') ; //falta dr

				//double eval_KN = es*eps*(1.0 - eta_obs*eta);
				double eval_KN = es*eps*(1.0 - mu);
				suma_2 += cLight*KNcrossSection(eval_KN)*nph*de;

				eps = eps*e_int;
			}

			double eval_e = es*Erest;
			double Ne = ne_eval(particle, eval_e, psc);
		
			suma_2 = suma_2*Ne / es;
	
			//--------------------------
			//continua la integral en eta
			double r = r_eta(eta, z, bG);
			double r2 = r_eta(eta - deta, z, bG);
			double dr = abs(r - r2);

			double x = sqrt(P2(r) + P2(z));
			double corr = P2(gamma*(x - bG*z))*x / (z*r);
			double corr2 = 2.0*r*h_d/P2(x);
			suma += (suma_1 + suma_2)*deta*(1.0 - mu)*dphi*corr*corr2; //dr

			phi = phi + dphi;
		}

		eta = eta + deta;
	}

	return suma*es*E; //ver por que va la segunda E
}



