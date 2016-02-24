#include "injection.h"

#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>

#include <fmath\physics.h>

//#include <fluminosities\luminositySynchrotron.h>
//#include <fluminosities\luminosityIC.h>


double primaryInjection(double E, double z, double t, Particle& particle);


void injection(Particle& p, State& st)
{
	double total = 0.0;
	p.ps.iterate([&p, &st, &total](const SpaceIterator& i){
		double factor = 0.0, E = i.par.E, z = i.par.R, t = i.par.T;
		switch (p.type)	{
		/*	case PT_photon:
				factor = radius / (volume*P2(E)*cLight); //este factor pasa de luminosidad a densidad de fotones [ 1/erg cm^3]

				total = 
					  factor*luminositySynchrotron(E, st.electron)
					+ factor*luminosityIC(E, st.electron, st.tpf);   //esto lo pongo para que se sumen todas las contribuciones
					
				total +=
					  factor*luminositySynchrotron(E, st.proton)
					+ factor*luminosityHadronic(E, st.proton)
					+ factor*luminosityPhotoHadronic(E, st.proton, st.tpf); 

			    break; */ 

			case PT_electron:
			//case PT_proton:
				total = primaryInjection(E, z, t, p);
				break;

		}
		p.injection.set(i, total);
	});
}


double nonThermalLuminosity(double z, double t)  //esta es la función que depende del número de estrellas a tiempo t
{
	return 1.0e43 / z;   //VER ingrese esto de prueba
}


double powerLaw(double E, double Emin, double Emax)        
{
	double result = pow(E, (-primaryIndex))*exp(-E / Emax)*exp(-5 * Emin / E);
	return result;
}
///////////////////synchr losses for secondary pairs

double normalization(double z, double t, double Emin, double Emax)
{

	double integral = RungeKuttaSimple(Emin, Emax, [&Emax, &Emin](double E){
		return E*powerLaw(E, Emin, Emax);
	});  //integra E*Q(E)  entre Emin y Emax

	double Q0 = nonThermalLuminosity(z, t) / integral;  //factor de normalizacion de la inyeccion
	return Q0;
}

double primaryInjection(double E, double z, double t, Particle& particle)
{
	ParticleType particleName = particle.type;

	double Emax = 1.6e-12*pow(10.0,particle.logEmax);
	double Emin = 1.6e-12*pow(10.0,particle.logEmin);

	double Q = powerLaw(E, Emin, Emax)*normalization(z, t, Emin, Emax);
	return Q;

}