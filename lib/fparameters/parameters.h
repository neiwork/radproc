#pragma once

#include <fmath\mathematics.h>

class Parameters {
public:
	double magneticField;
	double density;
	double radius;  //cambie Rc por radius, asi que no van a compilar las cosas viejas
	double volume;
	double accEfficiency;
	double starT;
	//double Rsp;

	//jet data
	double openingAngle;
	double B0;
	double Gamma;
	double Dlorentz;
	double Lj;

	double distance;

	double photonLogEmin;
	double photonLogEmax;

	//fun targetPhotonField;   //equivalente a double (*photonField)(double) es un puntero a funcion
	double targetPhotonEmin;
	double targetPhotonEmax;
	int nPhotonEnergies;
	

	// Data of electrons and protons

	double primaryIndex;
	double electronNormalization;
	double protonNormalization;

	double electronLogEmin;
	double electronLogEmax;
	double protonLogEmin;
	double protonLogEmax;

	double pionLogEmin;
	double pionLogEmax;
	double muonLogEmin;  //en realidad es log(2*mc^2) = 8.32
	double muonLogEmax;
	
	double neutrinoLogEmin;
	double neutrinoLogEmax;

	double neutronLogEmin;
	double neutronLogEmax;

	//Energy and time data

	double E;
	double R;
	double T;
	int nEnergies;        //massive particles
	int nTimes;

	double timeMin;
	double timeMax;   //segundos en dos horas


	double rmin;
	double rmax;
	int nR;
	
	Parameters();

};

extern Parameters parameters;