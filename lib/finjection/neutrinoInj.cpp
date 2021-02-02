#include "neutrinoInj.h"

//#include "dataInjection.h"
#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/interpolation.h>
#include <fmath/physics.h>
#include <algorithm>


//// las que usan muonNeuInj y muonAntiNeuInj son:
//
//double fQ_pion(double Ec, void* voiddata);          // esta es para inyectar neutinos muonicos a partir de los pi+
//double fQ_muon_L(double Ec, void* voiddata);   // esta es para inyectar neutinos muonicos a partir de los mu-L
//                                              //y antineutrinos muonicos a partir de los mu+L
//
//
//double fQ_muon_R(double Ec, void* voiddata);   // esta es para inyectar neutinos muonicos a partir de los mu-R
//                                              //y antineutrinos muonicos a partir de los mu+R
//
//// las que usan electronNeuInj y electronAntiNeuInj son:
//double fQ_elecNeu_L(double Ec, void* voiddata);  // esta es para inyectar neutinos electronicos a partir de los mu+L
//                                                 //y antineutrinos electronicos a partir de los mu-L
//
//double fQ_elecNeu_R(double Ec, void* voiddata); // esta es para inyectar neutinos electronicos a partir de los mu+R
//                                                 //y antineutrinos electronicos a partir de los mu-R
//
//
//

double Fmu_1(double x)
{
	return 5.0/3.0-3*x*x+4.0/3.0*x*x*x + (-1.0/3.0+3*x*x-8.0/3.0*x*x*x);
}

double Fmu_2(double x)
{
	return 5.0/3.0-3*x*x+4.0/3.0*x*x*x - (-1.0/3.0+3*x*x-8.0/3.0*x*x*x);
}

double Fmu_e1(double x)
{
	return 2.0-6*x*x+4*x*x*x + (2.0-12*x+18*x*x-8*x*x*x);
}

double Fmu_e2(double x)
{
	return 2.0-6*x*x+4*x*x*x - (2.0-12*x+18*x*x-8*x*x*x);
}

double muonNeutrinoInjection(double Enu, Particle& nu, Particle& mu, Particle& pi, const SpaceCoord& psc)
{
	double rpi = P2(muonMass/chargedPionMass);
	double inf = max(pi.emin(),Enu/(1.0-rpi));
	double injTau_Neutrino = integSimpsonLog(inf,pi.emax(),[Enu,rpi,&pi,&psc](double Epi)
								{
									double Npi = pi.distribution.interpolate({{0,Epi}},&psc);
									double tDecay = chargedPionMeanLife*(Epi/(chargedPionMass*cLight2));
									return Npi/tDecay/(Epi*(1.0-rpi));
								},30);
	inf = max(mu.emin(),Enu);
	double injMuon_Neutrino = 0.5*integSimpsonLog(inf,mu.emax(),[Enu,&mu,&psc](double Emu)
								{
									double Nmu = mu.distribution.interpolate({{0,Emu}},&psc);
									double tDecay = muonMeanLife*(Emu/(muonMass*cLight2));
									return Nmu/tDecay * (Fmu_1(Enu/Emu)+Fmu_2(Enu/Emu)) / Emu;
								},30);
	return injTau_Neutrino + injMuon_Neutrino;
}

double electronNeutrinoInjection(double Enu, Particle& nu, Particle& mu, const SpaceCoord& psc)
{
	double inf = max(mu.emin(),Enu);
	return 0.5*integSimpsonLog(inf,mu.emax(),[Enu,&mu,&psc](double Emu)
								{
									double Nmu = mu.distribution.interpolate({{0,Emu}},&psc);
									double tDecay = muonMeanLife*(Emu/(muonMass*cLight2));
									return Nmu/tDecay * (Fmu_e1(Enu/Emu)+Fmu_e2(Enu/Emu)) / Emu;
								},30);
}



//// las que usan muonNeuInj y muonAntiNeuInj son:

double h(double x)
{
	double r = P2(muonMass / chargedPionMass);

	return 2.0*r / ((1.0 - r)*x) - (1 + r) / (1 - r);
}

double fQ_pion(double E, double Ec, const Particle& creator, const SpaceCoord& psc)        // esta es para inyectar neutinos muonicos a partir de los pi+/-
{

	double distCreator = creator.distribution.interpolate({ { 0, Ec } }, &psc);
	//interpol(Ec,Ecreator,Ncreator,Ncreator.size()-1);  //el primer punto de Npi da feo

	double r = P2(muonMass / chargedPionMass);

	double x = E / Ec;

	double decayTime = chargedPionMeanLife*Ec / (chargedPionMass*cLight2);

	double Q_pion;

	//	if(1 > (r+x))	{  //esta condici�n es por la funci�n de Heaviside
	//esta comentada porque la transforme en una concidicion para el l�mite
	//inferior de la integrl
	Q_pion = distCreator / (Ec*(1 - r)*decayTime);
	//	}
	//	else { Q_pion = 0.0;	}

	return Q_pion;
}


double fQ_muon_L(double E, double Ec, const Particle& creator, const SpaceCoord& psc)   // esta es para inyectar neutinos muonicos a partir de los mu-L
																						//y antineutrinos muonicos a partir de los mu+L
{

	double distCreator = creator.distribution.interpolate({ { 0, Ec } }, &psc);
	//interpol(Ec,Ecreator,Ncreator,Ncreator.size()-1);

	double x = E / Ec;

	double decayTime = muonMeanLife*Ec / (muonMass*cLight2);

	double hache = -1.0; //h(x)

	double corchete = (5.0 / 3.0 - 3.0*P2(x) + 4.0*P3(x) / 3.0) + hache*(-1.0 / 3.0 + 3.0*P2(x) - 8.0*P3(x) / 3.0);

	double Q_muon = distCreator*corchete / (Ec*decayTime);

	return Q_muon;
}


double fQ_muon_R(double E, double Ec, const Particle& creator, const SpaceCoord& psc)   // esta es para inyectar neutinos muonicos a partir de los mu-R
																						//y antineutrinos muonicos a partir de los mu+R
{

	double distCreator = creator.distribution.interpolate({ { 0, Ec } }, &psc);
	//interpol(Ec,Ecreator,Ncreator,Ncreator.size()-1);

	double x = E / Ec;

	double decayTime = muonMeanLife*Ec / (muonMass*cLight2);

	double hache = 1.0; //h(x)

	double corchete = (5.0 / 3.0 - 3.0*P2(x) + 4.0*P3(x) / 3.0) + hache*(-1.0 / 3.0 + 3.0*P2(x) - 8.0*P3(x) / 3.0);

	double Q_muon = distCreator*corchete / (Ec*decayTime);

	return Q_muon;
}

double muonNeuInj(double E, int pos_E, int pos_t, const Particle& creator, const SpaceCoord& psc)
{
	std::string pName = creator.id;
	//ParticleType particleName = creator.type;

	double Emax = creator.emax();  
	/*st.photon.ps[DIM_R].size()
	int initial = pos_t*creator.ps[0].size();    //en estos dos debe ir proton y no particle
	int final   = (pos_t+1)*creator.ps[0].size();

	Vector Ncreator(creator.energyPoints.size(), 0.0);

	int j = 0;
	for (int i = initial; i < final; ++i)	{
		Ncreator[j] = creator.distribution.values[i];
		j = j+1;
	}*/

	double injection = 0.0;

	if (pName == "muon_L_minus") {
		//case PT_muon_L_minus:
		injection = RungeKuttaSimple(E, Emax, [&E, &Emax, &creator, &psc](double E) {
			return fQ_muon_L(E, Emax, creator, psc); });  //ver si Emax = Ec
		//	break;
	}
	else if (pName == "muon_R_minus") {
		//case PT_muon_R_minus:		
		injection = RungeKuttaSimple(E, Emax, [&E, &Emax, &creator, &psc](double E) {
			return fQ_muon_R(E, Emax, creator, psc); }); 
		//	break;
	}
	else if (pName == "pion") {
		//case PT_pion:
		double r = P2(muonMass / chargedPionMass);
		double Emin = std::max(E, E / (1.0 - r));
		injection = RungeKuttaSimple(Emin, Emax, [&E, &Emax, &creator, &psc](double E) {
			return fQ_pion(E, Emax, creator, psc); });  //multiplico por 0.5 porque supongo que tengo
																	   //mitad de pi+ -> muonNu y 
																	   //mitad de pi- -> antiMuonNu
	//	break;    
	}

	return injection;

}
//
//



//
//// las que usan electronNeuInj y electronAntiNeuInj son:
//double fQ_elecNeu_L(double Ec, void* voiddata);  // esta es para inyectar neutinos electronicos a partir de los mu+L
//                                                 //y antineutrinos electronicos a partir de los mu-L
//
//double fQ_elecNeu_R(double Ec, void* voiddata); // esta es para inyectar neutinos electronicos a partir de los mu+R
//                                                 //y antineutrinos electronicos a partir de los mu-R
//
//
//
//
//
//double muonAntiNeuInj(double E, int pos_E, int pos_t, Particle& creator)
//{
//	ParticleType particleName = creator.type;
//
//	double Emax = 1.6e-12*pow(10.0,creator.logEmax);  
//	
//	int initial = pos_t*creator.energyPoints.size();    //en estos dos debe ir proton y no particle
//	int final   = (pos_t+1)*creator.energyPoints.size();
//
//	Vector Ncreator(creator.energyPoints.size(), 0.0);
//	
//	int j = 0;
//	for (int i = initial; i < final; ++i)	{
//		Ncreator[j] = creator.distribution.values[i];
//		j = j+1;
//	}
//
//	DataInjection data;
//
//	data.E        = E;
//	data.Ncreator = Ncreator;
//	data.Ecreator = creator.energyPoints;
//
//	double injection = 0.0;
//
//	switch (particleName)	{
//		case PT_muon_L_plus: 
//			injection = RungeKuttaSimple(E, Emax, fQ_muon_L, &data);
//			break;
//		case PT_muon_R_plus:		
//			injection = RungeKuttaSimple(E, Emax, fQ_muon_R, &data);
//			break;
//		case PT_pion:
//			double r = P2(muonMass/chargedPionMass);
//			double Emin = std::max(E,E/(1.0-r));
//			injection = RungeKuttaSimple(Emin, Emax, fQ_pion, &data)*0.5;  //multiplico por 0.5 porque supongo que tengo
//			                                                               //mitad de pi+ -> muonNu y 
//																		   //mitad de pi- -> muonAntiNu
//			break;    
//	}
//
//	return injection;
//
//}
//
//
//
//
//
//double electronNeuInj(double E, int pos_E, int pos_t, Particle& creator)//el creator = muon
//{
//	ParticleType particleName = creator.type;
//
//	double Emax = 1.6e-12*pow(10.0,creator.logEmax);  
//	
//	int initial = pos_t*creator.energyPoints.size();    //en estos dos debe ir proton y no particle
//	int final   = (pos_t+1)*creator.energyPoints.size();
//
//	Vector Ncreator(creator.energyPoints.size(), 0.0);
//	
//	int j = 0;
//	for (int i = initial; i < final; ++i)	{
//		Ncreator[j] = creator.distribution.values[i];
//		j = j+1;
//	}
//
//	DataInjection data;
//
//	data.E        = E;
//	data.Ncreator = Ncreator;
//	data.Ecreator = creator.energyPoints;
//
//	double injection;
//
//	switch (particleName)	{
//		case PT_muon_L_plus: 
//			injection = RungeKuttaSimple(E, Emax, fQ_elecNeu_L, &data);
//			break;
//		case PT_muon_R_plus:		
//			injection = RungeKuttaSimple(E, Emax, fQ_elecNeu_R, &data);
//			break;
//	}
//		
//	return injection;
//
//}
//
//
//double electronAntiNeuInj(double E, int pos_E, int pos_t, Particle& creator)//el creator = muon
//{
//	ParticleType particleName = creator.type;
//
//	double Emax = 1.6e-12*pow(10.0,creator.logEmax);  
//	
//	int initial = pos_t*creator.energyPoints.size();    //en estos dos debe ir proton y no particle
//	int final   = (pos_t+1)*creator.energyPoints.size();
//
//	Vector Ncreator(creator.energyPoints.size(), 0.0);
//	
//	int j = 0;
//	for (int i = initial; i < final; ++i)	{
//		Ncreator[j] = creator.distribution.values[i];
//		j = j+1;
//	}
//
//	DataInjection data;
//
//	data.E        = E;
//	data.Ncreator = Ncreator;
//	data.Ecreator = creator.energyPoints;
//	
//	double injection;
//
//	switch (particleName)	{
//		case PT_muon_L_minus: 
//			injection = RungeKuttaSimple(E, Emax, fQ_elecNeu_L, &data);
//			break;
//		case PT_muon_R_minus:		
//			injection = RungeKuttaSimple(E, Emax, fQ_elecNeu_R, &data);
//			break;
//	}
//		
//	return injection;
//
//}
//
//
//
////definicion de funciones a integrar, las declaradas arriba



//
//
//
//
//
//
//double fQ_elecNeu_R(double Ec, void* voiddata)  // esta es para inyectar neutinos electronicos a partir de los mu+ 
//{
//	DataInjection* data = (DataInjection*)voiddata;
//	double E = data->E;      
//	Vector& Ncreator = data->Ncreator;
//	Vector& Ecreator = data->Ecreator;
//
//	double distCreator = interpol(Ec,Ecreator,Ncreator,Ncreator.size()-1);
//
//	double x = E/Ec;
//
//	double decayTime = muonMeanLife*Ec/(muonMass*cLight2);
//
//	double hache = 1.0; //h(x)
//
//	double corchete = (2.0-6*P2(x)+4*P3(x)) + hache*(2.0-12*x+18*P2(x)-8.0*P3(x));
//
//	double Q_muon = distCreator*corchete/(Ec*decayTime);
//
//	return Q_muon;
//}
//
//double fQ_elecNeu_L(double Ec, void* voiddata)    // esta es para inyectar antineutinos electronicos a partir de los mu- 
//{
//	DataInjection* data = (DataInjection*)voiddata;
//	double E = data->E;      
//	Vector& Ncreator = data->Ncreator;
//	Vector& Ecreator = data->Ecreator;
//
//	double distCreator = interpol(Ec,Ecreator,Ncreator,Ncreator.size()-1);
//
//	double x = E/Ec;
//
//	double decayTime = muonMeanLife*Ec/(muonMass*cLight2);
//
//	double hache = -1.0; //h(x)
//
//	double corchete = (2.0-6*P2(x)+4*P3(x)) + hache*(2.0-12*x+18*P2(x)-8.0*P3(x));
//
//	double Q_muon = distCreator*corchete/(Ec*decayTime);
//
//	return Q_muon;
//}
//
//
//
//
//
//
//
//
//
//
//
//
