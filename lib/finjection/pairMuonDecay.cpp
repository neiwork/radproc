#include "pairMuonDecay.h"

#include "dataInjection.h"
#include <fparameters/parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\interpolation.h>
#include <fmath\physics.h>
#include <algorithm>

double cMuonDec(double Gep, double E)
{
	double Ge = E;  //esto es para usar el paquete data

	double c  = Ge*Gep-sqrt((P2(Ge)-1)*(P2(Gep)-1));
		
	double Emin_mu = pow(10.0,muonLogEmin)*1.6e-12;

	double Gmu_min = Emin_mu/(muonMass*cLight2);
			
	double inf = std::max(c,Gmu_min);  //idem pero energia minima									   	

	return inf;
}

double dMuonDec(double Gep, double E)     //Gep=gama E prima
{
	double Ge = E;  //esto es para usar el paquete data

	double d = Ge*Gep+sqrt((P2(Ge)-1)*(P2(Gep)-1));

	double Emax_mu = pow(10.0,muonLogEmax)*1.6e-12;

	double Gmu_max = Emax_mu/(muonMass*cLight2);

	double sup = std::min(d,Gmu_max);  //esto no permite que el limite superior sea mayor que la energia max del muon
	
	return sup;               
}

double fMuonDec(double Gep, double Gmu, double E, Particle& creator)   //x=Ega; y=Eph
{ 

	double Tdec = muonMeanLife*E/(muonMass*cLight2);

	double Emu = Gmu*muonMass*cLight2;

	double Nmu = creator.dist(Emu);// interpol(Emu, Ecreator, Ncreator, Ncreator.size() - 1);

	double Qmu = Nmu*muonMass*cLight2/Tdec;  

	double Gep_max = 104;
		
	double P = 2*P2(Gep)*(3.0-2*Gep/Gep_max)/P3(Gep_max);

	double f = 0.5*P*Qmu/pow(((P2(Gep)-1.)*(P2(Gmu)-1.)),0.5);

	return f;

}



double pairMuonDecay(double E, Particle& particle, Particle& creator)    
{	
	using std::bind; using namespace std::placeholders; // para _1, _2, etc.
	

	double a = 1.0;
	double b = 104.0;

	double Erep	= electronMass*cLight2;

	E /= Erep;

	double integral  = RungeKutta(a,b,bind(cMuonDec,_1,E),bind(dMuonDec,_1,E),bind(fMuonDec,_1,_2,E,creator));

	



	double inj = integral/Erep;

	return inj;

}


	


	
	
