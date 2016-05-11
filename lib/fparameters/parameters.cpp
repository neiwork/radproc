#include "parameters.h"

Parameters::Parameters() {
	magneticField = 0.0;
	radius = 0.0;   //cambie Rc por radius, asi que no va a compilar las cosas viejas
	volume = 0.0;
	starT = 0.0;
	Lj = 0.0;

	//particle injection
	accEfficiency = 0.0;
	primaryIndex = 0.0;
	electronNormalization = 0.0;
	protonNormalization = 0.0;

	electronLogEmin = 0.0;
	electronLogEmax = 0.0;
	protonLogEmin = 0.0;
	protonLogEmax = 0.0;

	openingAngle = 0.0;
	Gamma		 = 0.0;
	Dlorentz	 = 0.0; //Dlorentz es el factor que transforma las distribuciones en el caso de jets

	photonLogEmin;
	photonLogEmax;
	nPhotonEnergies = 0;  //photons

	targetPhotonEmin = 0.0;
	targetPhotonEmax = 0.0;


	//Energy and time data

	nEnergies = 0;        //massive particles

	rmin = 0.0;
	rmax = 0.0;
	nR = 0;

	timeMin = 0.0;
	timeMax = 0.0;   //segundos en dos horas
	nTimes = 0;

	//particulas extras
	neutronLogEmin = 0.0;
	neutronLogEmax = 0.0;

	pionLogEmin = 0.0;
	pionLogEmax = 0.0;
	muonLogEmin = 0.0;  //en realidad es log(2*mc^2) = 8.32
	muonLogEmax = 0.0;

	neutrinoLogEmin = 0.0;
	neutrinoLogEmax = 0.0;
}

Parameters parameters;
