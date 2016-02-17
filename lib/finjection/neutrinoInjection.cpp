#include "neutrinoInjection.h"

#include "dataInjection.h"
#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\interpolation.h>
#include <fmath\physics.h>
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
//double muonNeuInj(double E, int pos_E, int pos_t, Particle& creator)
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
//		case PT_muon_L_minus:
//			injection = RungeKuttaSimple(E, Emax, fQ_muon_L, &data);
//			break;
//		case PT_muon_R_minus:		
//			injection = RungeKuttaSimple(E, Emax, fQ_muon_R, &data);
//			break;
//		case PT_pion:
//			double r = P2(muonMass/chargedPionMass);
//			double Emin = std::max(E,E/(1.0-r));  
//			injection = RungeKuttaSimple(Emin, Emax, fQ_pion, &data)*0.5;  //multiplico por 0.5 porque supongo que tengo
//			                                                               //mitad de pi+ -> muonNu y 
//																		   //mitad de pi- -> antiMuonNu
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
//double h(double x)
//{
//	double r = P2(muonMass/chargedPionMass);
//
//	return 2.0*r/((1.0-r)*x) - (1+r)/(1-r);
//}
//
//double fQ_pion(double Ec, void* voiddata)        // esta es para inyectar neutinos muonicos a partir de los pi+/-
//{
//	DataInjection* data = (DataInjection*)voiddata;
//	double E = data->E;      
//	Vector& Ncreator = data->Ncreator;
//	Vector& Ecreator = data->Ecreator;
//
//	double distCreator = interpol(Ec,Ecreator,Ncreator,Ncreator.size()-1);  //el primer punto de Npi da feo
//
//	double r = P2(muonMass/chargedPionMass);
//
//	double x = E/Ec;
//
//	double decayTime = chargedPionMeanLife*Ec/(chargedPionMass*cLight2);
//
//	double Q_pion;  
//
////	if(1 > (r+x))	{  //esta condición es por la función de Heaviside
//                       //esta comentada porque la transforme en una concidicion para el límite
//	                   //inferior de la integrl
//		Q_pion = distCreator/(Ec*(1-r)*decayTime);
////	}
////	else { Q_pion = 0.0;	}
//	
//	return Q_pion;   
//}
//
//
//double fQ_muon_L(double Ec, void* voiddata)   // esta es para inyectar neutinos muonicos a partir de los mu-L
//                                              //y antineutrinos muonicos a partir de los mu+L
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
//	double corchete = (5.0/3.0-3.0*P2(x)+4.0*P3(x)/3.0)+hache*(-1.0/3.0+3.0*P2(x)-8.0*P3(x)/3.0);
//
//	double Q_muon = distCreator*corchete/(Ec*decayTime);
//
//	return Q_muon;
//}
//
//
//double fQ_muon_R(double Ec, void* voiddata)   // esta es para inyectar neutinos muonicos a partir de los mu-R
//                                              //y antineutrinos muonicos a partir de los mu+R
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
//	double corchete = (5.0/3.0-3.0*P2(x)+4.0*P3(x)/3.0)+hache*(-1.0/3.0+3.0*P2(x)-8.0*P3(x)/3.0);
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
