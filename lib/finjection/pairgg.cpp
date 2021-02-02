#include "pairgg.h"


#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
//#include <algorithm>

double cAnnihilation(double x,double E)
{
	return x*P2(electronRestEnergy) / (4.0*E*(x-E));
}

double cAnnihilationPrueba(double x,double g)
{
	return x / (4.0*g*(x-g));
}

/*double dAnnihilation(double x )
{
	return targetPhotonEmax;        //este es el infinito del limite superior para la integral en Eph
}*/

double fAnnihilation(double Eg, double Eph, double Ee, const ParamSpaceValues& ntPh,
			const ParamSpaceValues& tpf, const SpaceCoord& distCoord, 
				double tpEminSoft, double tpEmaxSoft, double tpEminG, double tpEmaxG)   
{ 
	double Erest = electronRestEnergy;
	double photonDist_g(0.0), photonDist_soft(0.0);
	photonDist_g = (Eg > tpEminG && Eg < tpEmaxG) ?
		ntPh.interpolate({ { 0, Eg } }, &distCoord) : 0.0;
	photonDist_soft = (Eph > tpEminSoft && Eph < tpEmaxSoft) ?
		tpf.interpolate({ { 0, Eph } }, &distCoord) : 0.0;

	double result  =   photonDist_g*photonDist_soft / (Eg*Eg*Eg*Eph*Eph) *     
					   ( (4.0*Eg*Eg / (Ee*(Eg-Ee))) * log(4.0*Ee*Eph*(Eg-Ee)/(Erest*Erest*Eg))
					   - 8.0*Eg*Eph/(Erest*Erest) + 2.0*(2.0*Eg*Eph-Erest*Erest)*P2(Eg/Erest)/(Ee*(Eg-Ee))
					   - (1.0-Erest*Erest/(Eg*Eph)) * pow(Eg,4)/(P2(Ee*(Eg-Ee))) );
	return (Eph < Erest && Eph > P2(Erest)/Eg ) ? result : 0.0;   //pido que de algo solo si epsilon < Erep   //&& Erest<x
}


double fAnnihilationAux(double Eg, double Eph, double Ee)   
{ 
	double result = 0.0;
	double xg = Eg/electronRestEnergy;
	double xph = Eph/electronRestEnergy;
	double g = Ee/electronRestEnergy;
	if (xg > 1.0001) {
		result  = 4.0*xg*xg / (g*(xg-g)) * log(4.0*g*xph*(xg-g)/xg)
						   - 8.0*xg*xph + 2.0*(2.0*xg*xph-1.0)*xg*xg/(g*(xg-g))
						   - (1.0-1.0/(xg*xph)) * pow(xg,4)/P2(g*(xg-g));
	}
	return (xph*xg > 1.0) ? result : 0.0;   //pido que de algo solo si epsilon < Erep   //&& Erest<x
}

double fAnnihilationPrueba(double eg,double om,double g,const ParamSpaceValues& tpf,
							const SpaceCoord& distCoord,double tpEmin, double tpEmax)   
{ 
	double electronRestEnergy = electronMass*cLight2;
	double photonDist_x(0.0);
	if (om >= tpEmin/electronRestEnergy && om <= tpEmax/electronRestEnergy)
		photonDist_x = tpf.interpolate({{0,om*electronRestEnergy}},&distCoord);

	double result = photonDist_x/(om*om) *     
					   ( 4.0*eg*eg/(g*(eg-g)) * log(4.0*g*om*(eg-g)/eg)
					   - 8.0*eg*om + 2.0*eg*eg*(2.0*eg*om-1.0)
					   - (1.0-1.0/(eg*om))*eg*eg*eg*eg/(g*g*(eg-g)*(eg-g)) );
	return (om < 0.1) ? result : 0.0;
}


double pairGammaGamma(double Ee, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf,
					const SpaceCoord& distCoord, double tpEminSoft, double tpEmaxSoft,
					double tpEminG, double tpEmaxG)
{	
	using std::bind; using namespace std::placeholders; // para _1, _2, etc.
	double cte = 3.0*cLight*thomson*pow((electronMass*cLight2),4)/32.0;
	//double sup1 = tpEmaxG;  //este es el infinito del limite superior para la integral en Egamma
	//double inf = E;  //Ega_min < Ee_min  --> la condicion esta asegurada

	double integral  = RungeKutta(Ee,tpEmaxG, 
		[Ee,tpEminSoft](double Eg) {return max(tpEminSoft,cAnnihilation(Eg,Ee));},  //limite inferior
		[tpEmaxSoft](double Eg) {return tpEmaxSoft;},							  //limite superior
		[Ee,&ntPh,&tpf,&distCoord,tpEminSoft,tpEmaxSoft,tpEminG,tpEmaxG](double Eg, double Eph) 
		{return fAnnihilation(Eg,Eph,Ee,ntPh,tpf,distCoord,tpEminSoft,tpEmaxSoft,tpEminG,tpEmaxG); });
		
	double emissivityA = cte*integral;
	return emissivityA;
}

/*
void fggFunction(State& st, Vector& fVec, Vector& gVec, Vector& egVec, Vector& ephVec);
{
	FILE *fFile, *gFile, *egFile, *ephFile;
	fFile = fopen("f_gg.bin","rb");
	gFile = fopen("g.bin","rb");
	egFile = fopen("eg.bin","rb");
	ephFile = fopen("eph.bin","rb");
	
	float *fVec, *gVec, *egVec, *ephVec;
	fVec = (float*)calloc(100*100*100,sizeof(float));  
	gVec = (float*)calloc(100,sizeof(float)); 
	egVec = (float*)calloc(100,sizeof(float));
	ephVec = (float*)calloc(100,sizeof(float));
	fread(fVec,sizeof(float),100*100*100,fFile);
	fread(gVec,sizeof(float),100,gFile);
	fread(egVec,sizeof(float),100,egFile);
	fread(ephVec,sizeof(float),100,ephFile);
	fclose(fFile);
	fclose(gFile);
	fclose(egFile);
	fclose(ephFile);
}

double fNew(double Ee, double Eg, double Eph, Vector fVec, Vector gVec, Vector egVec, Vector ephVec)
{
	double g = Ee/electronRestEnergy;
	double eg = Eg/electronRestEnergy;
	double eph = Eph/electronRestEnergy;

	double logg = log10(g);
	if (logg > gVec[0] && logg < gVec[100-1]) {
		size_t pos_g = fbinarySearch(gVec,0,100,logg)+1;
		double logg1 = gVec[pos_g-1];
		double logg2 = gVec[pos_g];
		double logeg = log10(eg);
		if (logeg > egVec[0] && logeg < egVec[100-1]) {
			size_t pos_eg = fbinarySearch(egVec,0,100,logeg)+1;
			double logeg1 = egVec[pos_eg-1];
			double logeg2 = egVec[pos_eg];
			double logeph = log10(eph);
			if (logeph > ephVec[0] && logeph < ephVec[100-1]) {
				size_t pos_eph = fbinarySearch(ephVec,0,100,logeph)+1;
				double logeph1 = ephVec[pos_eph-1];
				double logeph2 = ephVec[pos_eph];
				
				double f111 = fVec[((pos_g-1)*100+(pos_eg-1))*100+(pos_eph-1)];
				double f112 = fVec[((pos_g-1)*100+(pos_eg-1))*100+pos_eph];
				double f121 = fVec[((pos_g-1)*100+pos_eg)*100+(pos_eph-1)];
				double f122 = fVec[((pos_g-1)*100+pos_eg)*100+pos_eph];
				double f211 = fVec[(pos_g*100+(pos_eg-1))*100+(pos_eph-1)];
				double f212 = fVec[(pos_g*100+(pos_eg-1))*100+pos_eph];
				double f221 = fVec[(pos_g*100+pos_eg)*100+(pos_eph-1)];
				double f222 = fVec[(pos_g*100+pos_eg)*100+pos_eph];
				
				double f11 = (f112-f111)/(logeph2-logeph1) * (logeph-logeph1) + f111;
				double f12 = (f122-f121)/(logeph2-logeph1)*(logeph-logeph1) + f121;
				double f1 = (f12-f11)/(logeg2-logeg1)*(logeg-logeg1) + f11;
				
				double f21 = (f212-f211)/(logeph2-logeph1) * (logeph-logeph1) + f211;
				double f22 = (f222-f221)/(logeph2-logeph1) * (logeph-logeph1) + f221;
				double f2 = (f22-f21)/(logeg2-logeg1)*(logeg-logeg1) + f21;
				
				double fa = (f2-f1)/(logg2-logg1)*(logg-logg1) + f1;
				return pow(10.0,fa);
			} else
				return 0.0;
		} else
			return 0.0;
	} else
		return 0.0;
}
*/
double pairGammaGammaNew(double Ee, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf,
					const SpaceCoord& distCoord, double tpEminSoft, double tpEmaxSoft,
					double tpEminG, double tpEmaxG)
{	
	double cte = 3.0*cLight*thomson*pow((electronMass*cLight2),4)/32.0;
	double integral  = integSimpsonLog(Ee*1.1,tpEmaxG,
			[Ee,&ntPh,&tpf,&distCoord,tpEminSoft,tpEmaxSoft,tpEminG,tpEmaxG](double Eg)
			{
				double inf = max(tpEminSoft,cAnnihilation(Eg,Ee));
				double ntph = (Eg > tpEminG && Eg < tpEmaxG) ? ntPh.interpolate({{0,Eg}},&distCoord) : 0.0;
				double integ1 = integSimpsonLog(inf,tpEmaxSoft,
						[Eg,Ee,&tpf,&distCoord,tpEminSoft,tpEmaxSoft](double Eph)
						{
							double nPh = (Eph > tpEminSoft && Eph < tpEmaxSoft) ?
									tpf.interpolate({{0,Eph}},&distCoord) : 0.0;
							return fAnnihilationAux(Eg,Eph,Ee)*nPh/(Eph*Eph);
						},50);
				return ntph*integ1/(Eg*Eg*Eg);
			},50);
	double emissivityA = cte*integral;
	return emissivityA;
}

/*
double pairGammaGammaInterp(double Ee, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf,
					const SpaceCoord& distCoord, double tpEminSoft, double tpEmaxSoft,
					double tpEminG, double tpEmaxG)
{	
	double cte = 3.0*cLight*thomson*pow((electronMass*cLight2),4)/32.0;
	double integral  = integSimpsonLog(Ee*1.1,tpEmaxG,
			[Ee,&ntPh,&tpf,&distCoord,tpEminSoft,tpEmaxSoft,tpEminG,tpEmaxG](double Eg)
			{
				double inf = max(tpEminSoft,cAnnihilation(Eg,Ee));
				double ntph = (Eg > tpEminG && Eg < tpEmaxG) ? ntPh.interpolate({{0,Eg}},&distCoord) : 0.0;
				double integ1 = integSimpsonLog(inf,tpEmaxSoft,
						[Eg,Ee,&tpf,&distCoord,tpEminSoft,tpEmaxSoft](double Eph)
						{
							double nPh = (Eph > tpEminSoft && Eph < tpEmaxSoft) ?
									tpf.interpolate({{0,Eph}},&distCoord) : 0.0;
							return fNew(Ee,Eg,Eph)*nPh/(Eph*Eph);
						},50);
				return ntph*integ1/(Eg*Eg*Eg);
			},50);
	double emissivityA = cte*integral;
	return emissivityA;
}
*/

double pairInjectionPrueba(double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf,
					const SpaceCoord& distCoord, double tpEmin, double tpEmax)
{	
	double electronRestEnergy = electronMass*cLight2;
	double gamma = E / electronRestEnergy;
	double cte = 3.0*cLight*thomson*electronRestEnergy/32.0;
	double sup = tpEmax/electronRestEnergy;  //este es el infinito del limite superior para la integral en Egamma

	double integral  = RungeKuttaSimple(gamma,sup,[&](double eg){ return 
						RungeKuttaSimple(cAnnihilationPrueba(eg,gamma),tpEmax/electronRestEnergy,[&](double om) { return
						((eg >= tpEmin/electronRestEnergy && eg <= tpEmax/electronRestEnergy) ?
						ntPh.interpolate({{0,eg*electronRestEnergy}},&distCoord) : 0.0) / (eg*eg*eg) *
						fAnnihilationPrueba(eg,om,gamma,tpf,distCoord,tpEmin,tpEmax);});});
		
	double emissivityA = cte*integral;
	return 2.0*emissivityA;
}