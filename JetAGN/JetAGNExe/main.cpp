#include <JetAGN/jetAGN.h>

int main() {
	return jetAGN();
}

/*	int useCase(){

	coord system = cilinder
	physical dimensions = 1

	injected particles:
	  electron:
	observed particles:
      photon:
	radprocs:
	  synchrotron ( particle=electron, field=magfield )=>photon
	  inverse compton ( particla=electron, field=photons(target) )=>photons


	----------

	inject electron (electron.injection = ...)

	distribution(list(procs))

	distribution
		electron.distribution = distribution(injection) {
			* evoluciona esa inyeccion inicial a un estado estable
			* radprocs hacen contribuciones via terminos en una eq diferencial a resolver:
			  synchrotron:
			    electron.distribution -= synchrotronLosses(...)
				ic.distribution -= icLosses(...)

	{ distribucion de electrones esta lista }

	photon.distribution +=
		contribuciones radproc:
		  synchrotronEmi()
		  icEmi()

	emissivity(procs)

	-------------------------


	radproc:
	  map<PSV> => map<PSV>
	 






	 } */