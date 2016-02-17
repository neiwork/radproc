#include "parameters.h"

Parameters::Parameters() {
	magneticField = 0.0;
	density = 0.0;
	radius = 0.0;   //cambie Rc por radius, asi que no va a compilar las cosas viejas
	volume = 0.0;


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
	B0			 = 0.0;
	Gamma		 = 0.0;
	Dlorentz	 = 0.0; //Dlorentz es el factor que transforma las distribuciones en el caso de jets

	distance = 0.0;

	photonLogEmin;
	photonLogEmax;
	nPhotonEnergies = 0;  //photons


	//fun targetPhotonField   = NULL;   //equivalente a (*photonField)(double) es un puntero a funcion
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


double magneticField = 0.0;
double density = 0.0;
double radius = 0.0;   //cambie Rc por radius, asi que no va a compilar las cosas viejas
double volume = 0.0;
double accEfficiency = 0.0;

double openingAngle = 0.0;
double B0 = 0.0;
double Gamma = 0.0;
double Dlorentz = 0.0; //Dlorentz es el factor que transforma las distribuciones en el caso de jets

double distance = 0.0;

double photonLogEmin;
double photonLogEmax;


//fun targetPhotonField   = NULL;   //equivalente a double (*photonField)(double) es un puntero a funcion
double targetPhotonEmin = 0.0;
double targetPhotonEmax = 0.0;

double primaryIndex = 0.0;
double factor_qrel = 0.0;
double steadyElectronNormalization = 0.0;
double steadyProtonNormalization = 0.0;   // es decir, a Lrel = 1% Lc
double electronNormalization = 0.0;
double protonNormalization = 0.0;

double electronLogEmin = 0.0;
double electronLogEmax = 0.0;
double protonLogEmin = 0.0;
double protonLogEmax = 0.0;



//Energy and time data

int nEnergies = 0;        //massive particles
int nPhotonEnergies = 0;  //259;  //photons
int nTimes = 0;

double rmin = 0.0;
double rmax = 0.0;
int nR = 0;

double timeMin = 0.0;
double timeMax = 0.0;   //segundos en dos horas
double timeStep = 0.0;

//particulas extras
double neutronLogEmin = 0.0;
double neutronLogEmax = 0.0;

double pionLogEmin = 0.0;
double pionLogEmax = 0.0;
double muonLogEmin = 0.0;  //en realidad es log(2*mc^2) = 8.32
double muonLogEmax = 0.0;

double neutrinoLogEmin = 0.0;
double neutrinoLogEmax = 0.0;
