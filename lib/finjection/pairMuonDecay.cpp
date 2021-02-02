#include "pairMuonDecay.h"


#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
//#include <fmath\interpolation.h>
#include <fmath/physics.h>
//#include <algorithm>

double cMuonDec(double Gep, double Ge)
{
	//double Ge = E/(electronmass*cLight2); 

	double Gmu_min  = Ge*Gep-sqrt((P2(Ge)-1.0)*(P2(Gep)-1.0));

	return Gmu_min;
}

double dMuonDec(double Gep, double Ge)     //Gep=gama E prima
{
	///double Ge = E/(electronmass*cLight2); 

	double Gmu_max = Ge*Gep+sqrt((P2(Ge)-1.0)*(P2(Gep)-1.0));

	return Gmu_max;               
}

double fMuonDec(double Gep, double Gmu, Particle& c, const SpaceCoord& distCoord)   //x=Gep; y=Gmu; c = muon
{ 	
	
	double Emu = Gmu*muonMass*cLight2;
	
	double Tdec = muonMeanLife*Emu/(muonMass*cLight2);
	
	double Nmu = 0.0;

	if (Emu > c.emin() && Emu < c.emax()){
		Nmu = c.distribution.interpolate({ { 0, Emu} }, &distCoord);
	}
	
	double Qmu = Nmu*muonMass*cLight2/Tdec;  //N(G) = N(E)*mc^2

	double Gep_max = 104;
		
	double P = 2.0*P2(Gep)*(3.0-2.0*Gep/Gep_max)/P3(Gep_max);
	double f = 0.5*P*Qmu/sqrt((P2(Gep)-1.0)*(P2(Gmu)-1.0));

	return f;

}

double pairMuonDecay(double E, Particle& c, const SpaceCoord& distCoord)
{	
	double inf = 2.0;
	double sup = 104.0;

	double Gmu_min = c.emin()/(muonMass*cLight2);
	double Gmu_max = c.emax()/(muonMass*cLight2);
	
	double Ge = E /(electronMass*cLight2); //aparece solo en los limites

		double integral  = RungeKutta(inf,sup, 
		[Ge,Gmu_min](double Gep) {return max(Gmu_min,cMuonDec(Gep, Ge)); },  //limite inferior
		[Ge,Gmu_max](double Gep) {return max(Gmu_max,dMuonDec(Gep, Ge)); },	 //limite superior
		[&c,&distCoord](double Gep, double Gmu) {return fMuonDec(Gep, Gmu, c, distCoord); });
		
	//double integral  = RungeKutta(a,b,bind(cMuonDec,_1,E),bind(dMuonDec,_1,E),bind(fMuonDec,_1,_2,E,creator));


	double inj = integral/(electronMass*cLight2);

	return inj;

}

double pairMuonDecayNew(double Eee, Particle& c, const SpaceCoord& distCoord)
{	
	double inf = 1.0;
	double sup = 104.0;
	double gg = Eee/electronRestEnergy;
	double ss = sqrt(gg*gg-1.0);
	
	double integ1 = integSimpson(0.0,1.0,[gg,ss,sup,&c,&distCoord](double sp)
		{
			double gg_p = sqrt(sp*sp+1.0);
			double Pfun = 2.0*gg_p*gg_p/P3(sup) * (3.0-2.0*gg_p/sup);
			double inf2 = gg*gg_p-ss*sp;
			double sup2 = gg*gg_p+ss*sp;
			double infgm = c.emin()/(muonMass*cLight2);
			double supgm = c.emax()/(muonMass*cLight2);
			double integ11 = 0.0;
			double integ12 = 0.0;
			double spMax = sqrt(P2(min(sup2,supgm))-1.0);
			if (inf2 < infgm) {
				integ11 = integSimpson(0.0,1.0,[&c,&distCoord](double t)
				{
					double gmu = sqrt(t*t+1.0);
					double Emu = gmu*muonMass*cLight2;
					double Nmu = (Emu < c.emax()) ? c.distribution.interpolate({{0,Emu}},&distCoord) : 0.0;
					double tdecay = muonMeanLife*gmu;
					double Qmu = Nmu / tdecay;
					return Qmu / gmu;
				},5);
				integ12 = integSimpson(0.0,log(spMax),[&c,&distCoord](double logt)
				{
					double gmu = sqrt(exp(2.0*logt)+1.0);
					double Emu = gmu*muonMass*cLight2;
					double Nmu = (Emu < c.emax()) ? c.distribution.interpolate({{0,Emu}},&distCoord) : 0.0;
					double tdecay = muonMeanLife*gmu;
					double Qmu = Nmu / tdecay;
					return exp(logt) * Qmu / gmu;
				},30);
			} else {
				double spMin = sqrt(P2(inf2)-1.0);
				integ12 = integSimpson(log(spMin),log(spMax),[&c,&distCoord](double logt)
				{
					double gmu = sqrt(exp(2.0*logt)+1.0);
					double Emu = gmu*muonMass*cLight2;
					double Nmu = (Emu < c.emax()) ? c.distribution.interpolate({{0,Emu}},&distCoord) : 0.0;
					double tdecay = muonMeanLife*gmu;
					double Qmu = Nmu / tdecay;
					return exp(logt) * Qmu / gmu;
				},30);
			}
			return (integ11+integ12) * Pfun / gg_p;
		},5);
	double spMax = sqrt(P2(sup)-1.0);
	double integ2 = integSimpson(0.0,log(spMax),[gg,ss,sup,&c,&distCoord](double logsp)
		{
			double sp = exp(logsp);
			double gg_p = sqrt(sp*sp+1.0);
			double Pfun = 2.0*gg_p*gg_p/P3(sup) * (3.0-2.0*gg_p/sup);
			double inf2 = gg*gg_p-ss*sp;
			double sup2 = gg*gg_p+ss*sp;
			double infgm = c.emin()/(muonMass*cLight2);
			double supgm = c.emax()/(muonMass*cLight2);
			double integ11 = 0.0;
			double integ12 = 0.0;
			double spMax = sqrt(P2(min(sup2,supgm))-1.0);
			if (inf2 < infgm) {
				integ11 = integSimpson(0.0,1.0,[&c,&distCoord](double t)
				{
					double gmu = sqrt(t*t+1.0);
					double Emu = gmu*muonMass*cLight2;
					double Nmu = (Emu < c.emax()) ? c.distribution.interpolate({{0,Emu}},&distCoord) : 0.0;
					double tdecay = muonMeanLife*gmu;
					double Qmu = Nmu / tdecay;
					return Qmu / gmu;
				},5);
				integ12 = integSimpson(0.0,log(spMax),[&c,&distCoord](double logt)
				{
					double gmu = sqrt(exp(2.0*logt)+1.0);
					double Emu = gmu*muonMass*cLight2;
					double Nmu = (Emu < c.emax()) ? c.distribution.interpolate({{0,Emu}},&distCoord) : 0.0;
					double tdecay = muonMeanLife*gmu;
					double Qmu = Nmu / tdecay;
					return exp(logt) * Qmu / gmu;
				},30);
			} else {
				double spMin = sqrt(P2(inf2)-1.0);
				integ12 = integSimpson(log(spMin),log(spMax),[&c,&distCoord](double logt)
				{
					double gmu = sqrt(exp(2.0*logt)+1.0);
					double Emu = gmu*muonMass*cLight2;
					double Nmu = (Emu < c.emax()) ? c.distribution.interpolate({{0,Emu}},&distCoord) : 0.0;
					double tdecay = muonMeanLife*gmu;
					double Qmu = Nmu / tdecay;
					return exp(logt) * Qmu / gmu;
				},30);
			}
			return sp * (integ11+integ12) * Pfun / gg_p;
		},30);
	return (muonMass/electronMass)*(integ1+integ2);
}

double pairMuonDecayNew2(double Eee, Particle& c, const SpaceCoord& distCoord)
{	
	double factor = 3.0;
	double Qmu = (Eee*factor > c.emin() && Eee*factor < c.emax()) ?
					c.injection.interpolate({{0,Eee*factor}},&distCoord) : 0.0;
	return factor*2.0*Qmu;
}

