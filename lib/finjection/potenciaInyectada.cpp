#include "potenciaInyectada.h"


#include <fluminosities\factorConversion.h>
#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\interpolation.h>


double fPotIny(double Ee, Particle& c)
{
	Vector& Einy = c.ps.dimensions[0]->values;
	Vector& Qiny = c.injection.values;

	//DataInjection* data = (DataInjection*)voiddata;
	//Vector& Qiny = data->Ncreator;
	//Vector& Einy = data->Ecreator;

	double iny = interpol(Ee, Einy, Qiny, Qiny.size() - 1);

	return Ee*iny;
}

double fPotInyPhoton(double Ee, Particle& c)
{
	Vector& Einy = c.ps.dimensions[0]->values;
	Vector& Qiny = c.injection.values;

	//DataInjection* data = (DataInjection*)voiddata;
	//Vector& Qiny = data->Ncreator;
	//Vector& Einy = data->Ecreator;

	double iny = interpol(Ee, Einy, Qiny, Qiny.size() - 1);

	double factor = factorLumToQph(Ee); //este factor pasa de luminosidad a
	//densidad de fotones [ 1/erg cm^3] 
	double lum = iny*factor;

	return iny / Ee;
}


double potenciaInyectada(Particle& p) 
{                  
	ParticleType particleName = p.type; 

	double Emax = 1.6e-12*pow(10.0,p.logEmax);  
	double Emin = 1.6e-12*pow(10.0,p.logEmin);  

	//DataInjection data;
	//data.Ncreator  = particle.injection.values;
	//data.Ecreator  = particle.energyPoints;

	fun1 rk_fQl = [&p](double e){
		return fPotIny(e, p);
	};

	double integral = RungeKuttaSimple(Emin, Emax, rk_fQl);
	//double integral = RungeKuttaSimple(Emin, Emax, fPotIny, &data); 

	double injection = integral*volume;
	
	return injection;

}


//	if (particleName == PT_photon){
//		integral = RungeKuttaSimple(Emin, Emax, fPotInyPhoton, &data); 
//		injection = integral*volume;
//	}
//	else{
//		integral = RungeKuttaSimple(Emin, Emax, fPotIny, &data); 	
//	}
	



//double potenciaInyectada2(Vector& energy, Vector& injection) 
//{
//	double Emax = 1.6e-12*targetPhotonEmax;  
//	double Emin = 1.6e-12*targetPhotonEmin;  
//
//	DataInjection data;
//	data.Ncreator  = injection;
//	data.Ecreator  = energy;
//
//	double integral = RungeKuttaSimple(Emin, Emax, fPotIny, &data); 
//
//	double inj = integral*volume;
//	
//	return inj;
//
//}
