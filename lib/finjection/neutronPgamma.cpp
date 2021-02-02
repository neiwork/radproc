#include "neutronPgamma.h"


#include "pgammaPionInj.h"

#include <fparameters/parameters.h>

//con incluir #include "pgammaPionInj.h" ya incluyo las siguientes funciones

//double t_pion_PH(double E, Particle& particle, fun1 tpf);   
//double omegaPH(double E, Particle& particle, fun1 tpf);



double psiEsc(double E, Particle& neutron, const ParamSpaceValues& tpf, const SpaceCoord& distCoord,
				double tpEmin, double tpEmax)  //E=Ep
{

	double seda_pn = 0.5;  //probabilidad de que el n se convierta en proton psi_np = 1-psi_nn

	double tauDec = neutronMeanLife*E/(neutronMass*cLight2);

	double L = 0.0;

	int nT = 10;
	
	double tcross = 0.0;
	double dt = tcross/nT;
	double t = 0.0;

	for (size_t i=0; i < nT; ++i){

		L += tcross*(seda_pn*omegaPH(E,neutron,tpf,distCoord,tpEmin,tpEmax)+1.0/tauDec)*dt;

		t += dt;

	}

	double result = exp(-L);

	return result; //psiEsc*psi_pn*omegaPH(E, particle, tpf, distCoord, tpEmin, tpEmax);
}


double neutronPgamma(double En, Particle& neutron, Particle& proton, const ParamSpaceValues& tpf,
						const SpaceCoord& distCoord, double tpEmin, double tpEmax)  
{
	//double protonDist = proton.dist(E);// interpol(E, proton.energyPoints, Nproton, Nproton.size() - 1);
	//double t_1 = t_pion_PH(E,proton,tpf,distCoord,tpEmin,tpEmax);     //esto no es lossesPH porque son perdidas solo del canal de produccion de piones
	//double omega_pg = omegaPH(E,proton,tpf,distCoord,tpEmin,tpEmax);
	//double averageInel = t_1/omega;
	double averageInel = 0.3;
	double Ep = En/(1.0-averageInel);

	double xi_pn = 0.5;
	double protonDist = (Ep > proton.emin() && Ep < proton.emax()) ?
						proton.distribution.interpolate({{0,Ep}},&distCoord) : 0.0;

//		double k2 = 0.6;
//		double p1 = (k2-averageInel)/(k2-k1);	
//		double nChargedPion = 2.0-1.5*p1;

		//double pEsc = psiEsc(Ep, tcross, neutron, tpf, distCoord, tpEmin, tpEmax);
		//double tauEsc = 1.0/(pEsc*seda_pn*omega);
	double omega_pg = omegaPH(Ep,proton,tpf,distCoord,tpEmin,tpEmax);
	double emissivity = xi_pn * protonDist / (1.0-averageInel) * omega_pg; // /tauEsc;
	return emissivity;
}