#include <stdio.h>

#include "state.h"
#include "modelParameters.h"
#include "radiativeLosses.h"
#include "injection.h"
////#include "photonInjection.h"
#include "distribution.h"
//#include "processes.h"
////#include "photonDistribution.h"

#include "write.h"
//#include "messages.h"
#include "checks.h"

//#include <fmath\matrixInit.h>

#include <fparticle\Particle.h>
#include <fparameters\parameters.h>
#include <fmath\physics.h>

//void injectionAndDistribution(Particle& p, State& st)
//{
//	injection(p, st);
//	distribution(p, st);
//}

int main ()
{
	//testBinarySearch();


	setParameters();

	//show_message(msgStart, Module_Main);

	State model;

	model.photon.injection.ps.iterate([&model](const SpaceIterator& i){
		double nph = blackBody(i.par.E, i.par.R);
		model.nph.set(i, nph);
	});


//	radiativeLosses(model);

	injection(model.electron, model);
//	writeAllSpaceParam("electronInj.txt", model.electron.injection);
	writeEnergyFunction("electronInj_E.txt", model.electron.injection, 1, 1); //escribe Q(E), para r(0) y t(0)
	
	distribution(model.electron, model);
	writeAllSpaceParam("electronDist.txt", model.electron.distribution);
	writeEnergyFunction("electronDist_E.txt", model.electron.distribution, 1, 1);


	//injection(model.photon, model);  //aca calculo las luminosidades

	return 0;
}






	//show_message(msgStart, Module_primaryInjection);

	////injectionAndDistribution(model.electron, model);
	//injectionAndDistribution(model.proton, model);
	//writeDistribution("protonDist.txt", model.proton);

	//show_message(msgEnd, Module_primaryInjection);

	//// Luminosidad primaria
	//show_message(msgStart, Module_primaryLuminosity); 
	//injection(model.photon, model);
	//show_message(msgEnd, Module_primaryLuminosity); 
	//
	//// Distribucion fotones
	////photondistribution(model);

	////writeLum("finalLuminosityDis.txt", st.photon.injection, st.photon);
	////writeFlux("finalFlux.txt", st.photon.injection, st.photon);

	//processes(model);

	//show_message(msgEnd, Module_Main);




