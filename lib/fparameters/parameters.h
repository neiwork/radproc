#pragma once

#include <fmath\mathematics.h>

extern double magneticField;
extern double density;
extern double radius;  //cambie Rc por radius, asi que no van a compilar las cosas viejas
extern double volume;
extern double accEfficiency;

extern double openingAngle;
extern double B0;
extern double Gamma;
extern double Dlorentz;

extern double distance;
//extern double tcross;

extern double photonLogEmin;
extern double photonLogEmax;

//extern fun targetPhotonField;   //equivalente a double (*photonField)(double) es un puntero a funcion
extern double targetPhotonEmin;
extern double targetPhotonEmax;
extern int nPhotonEnergies;



// Data of electrons and protons

extern double primaryIndex;
extern double electronNormalization;
extern double protonNormalization;

extern double electronLogEmin;
extern double electronLogEmax;
extern double protonLogEmin;
extern double protonLogEmax;


//Energy, zPoints and time data

extern int nEnergies;        //massive particles

extern double rmin;
extern double rmax;
extern int nR;

extern double timeMin;
extern double timeMax;   //segundos en dos horas
extern int nTimes;


//particulas extras
extern double neutronLogEmin;
extern double neutronLogEmax;

extern double pionLogEmin;
extern double pionLogEmax;
extern double muonLogEmin;  //en realidad es log(2*mc^2) = 8.32
extern double muonLogEmax;

extern double neutrinoLogEmin;
extern double neutrinoLogEmax;

extern double neutronLogEmin;
extern double neutronLogEmax;



class Parameters {
public:
	double magneticField;
	double density;
	double radius;  //cambie Rc por radius, asi que no van a compilar las cosas viejas
	double volume;
	double accEfficiency;

	//jet data
	double openingAngle;
	double B0;
	double Gamma;
	double Dlorentz;

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