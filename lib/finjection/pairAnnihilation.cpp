#include "pairAnnihilation.h"

#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\interpolation.h>
#include <fmath\physics.h>

double gammaU(double k, double gamma_min, double gamma_mas)
{
	if( (k > gamma_mas && k < gamma_min) || (k < gamma_mas && k > gamma_min) )
	{
		double beta_mas = sqrt(1.0-1/P2(gamma_mas));
		double beta_min = sqrt(1.0-1/P2(gamma_min));
		
		return sqrt(0.5*(1.0+gamma_mas*gamma_min*(1+beta_mas*beta_min)));
	}
	else
	{
		return k*(gamma_mas+gamma_min-k);
	}
}

double cPairAnn(double x, double E)
{	
	//AnnihilationData* data = (AnnihilationData*)voiddata;
	//const double E = data->E;

	return P2(E)/x;   //1.6e-12*pow(10.0,electronLogEmin);
	//return P2(electronMass*cLight2)/x;
	
}

double dPairAnn(double x )
{
	return 1.6e-12*pow(10.0,protonLogEmax);               
}

double fPairAnn(double x, double y, double E, Particle& priElectron, Particle& secElectron, Particle& positron)   //x=Ee; y=Ep positron
{ 
	//AnnihilationData* data = (AnnihilationData*)voiddata;
	//const double E = data->E;
	//const Vector& Ee(data->Ee);   //secondaryElectrons
	//const Vector& Ne(data->Ne);
	//const Vector& Eep(data->Eep);   //primaryElectrons
	//const Vector& Nep(data->Nep);
	//const Vector& Ep(data->Ep);   //positrons
	//const Vector& Np(data->Np);

	double distElectrons = priElectron.dist(x) + secElectron.dist(x);// (x, Ee, Ne, Ne.size() - 1) + interpol(x, Eep, Nep, Nep.size() - 1);
	double distPositrons = positron.dist(y);// interpol(y, Ep, Np, Ne.size() - 1);

	//double aux = x*y/P2(electronMass*cLight2);

	double k = E/(electronMass*cLight2);
	double gamma_min = x/(electronMass*cLight2);
	double gamma_mas = y/(electronMass*cLight2);

	double parenthesis =  P2(gammaU(k,gamma_min,gamma_mas))/(std::abs(k-gamma_mas)+2.0/pi)
		                + P2(gammaU(k,gamma_min,gamma_mas))/(std::abs(k-gamma_min)+2.0/pi);

	double result = distElectrons*distPositrons*parenthesis/P2(x*y);

//	return x*y >= P2(E)  ? result : 0.0;  //lo paso a limite de integracion
	return result;

}

double pairAnnihilation(double E, Particle& electron, Particle& secondaryElectron, Particle& positron)    
{	
	using std::bind; using namespace std::placeholders; // para _1, _2, etc.

	double cte = 3 * cLight*thomson*P3(electronMass*cLight2) / 8.0;

	//double supEmax = 1.6e-12*pow(10.0,secondaryElectron.logEmax);  
	//double infEmin = 1.6e-12*pow(10.0,secondaryElectron.logEmin);
	//
	//AnnihilationData adata;    //constructor

	//adata.E    = E;
	//adata.Eep  = electron.energyPoints;
	//adata.Nep = electron.distribution.values;
	//adata.Ee   = secondaryElectron.energyPoints;
	//adata.Ne = secondaryElectron.distribution.values;
	//adata.Ep   = positron.energyPoints;
	//adata.Np = positron.distribution.values;

	double integral = RungeKutta(
		secondaryElectron.emin(), 
		secondaryElectron.emax(), 
		bind(cPairAnn, _1, E), dPairAnn,
		bind(fPairAnn, _1, _2, E, electron, secondaryElectron, positron));

	double emissivity = cte*integral;

	return emissivity;

}

//aca empieza lo que tenia de coppi y Blandford, pero con esa integral el indice siempre es Eph^-1

//double cPairAnn(double x, void* voiddata)
//{	
//	return 1.6e-12*pow(10.0,electronLogEmin);
//	//return P2(electronMass*cLight2)/x;
//	
//}
//
//double dPairAnn(double x, void* )
//{
//	return 1.6e-12*pow(10.0,protonLogEmax);               
//}
//
//double fPairAnn(double x, double y, void* voiddata)   //x=Ee; y=Ep positron
//{ 
//	AnnihilationData* data = (AnnihilationData*)voiddata;
//	const double E = data->E;
//	const Vector& Ee(data->Ee);   //secondaryElectrons
//	const Vector& Ne(data->Ne);
//	const Vector& Eep(data->Eep);   //primaryElectrons
//	const Vector& Nep(data->Nep);
//	const Vector& Ep(data->Ep);   //positrons
//	const Vector& Np(data->Np);
//
//	double distElectrons = interpol(x,Ee,Ne,Ne.size()-1)+interpol(x,Eep,Nep,Nep.size()-1);
//	double distPositrons = interpol(y,Ep,Np,Ne.size()-1);
//
//	double aux = x*y/P2(electronMass*cLight2);
//
//	double result = (x+y)*(log(aux)+pow(aux,-0.5))*distElectrons*distPositrons/(x*y);
//
//	return x*y >= P2(E)  ? result : 0.0;   
//
//}
//
//double pairAnnihilation(double E, Particle& electron, Particle& secondaryElectron, Particle& positron)    
//{	
//	double cte = 3*cLight*thomson*P2(electronMass*cLight2)/8.0;
//
//	double supEmax = 1.6e-12*pow(10.0,secondaryElectron.logEmax);  
//	double infEmin = 1.6e-12*pow(10.0,secondaryElectron.logEmin);
//	
//	AnnihilationData adata;    //constructor
//
//	adata.E    = E;
//	adata.Eep  = electron.energyPoints;
//	adata.Nep  = electron.distribution;
//	adata.Ee   = secondaryElectron.energyPoints;
//	adata.Ne   = secondaryElectron.distribution;
//	adata.Ep   = positron.energyPoints;
//	adata.Np   = positron.distribution;
//
//	double integral  = RungeKutta(infEmin,supEmax,cPairAnn,dPairAnn,fPairAnn,&adata);
//
//	double emissivity = cte*integral/P2(E);
//
//	return emissivity;
//
//}