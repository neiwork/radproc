#pragma once // evita que el preprocesador incluya 
             // varias veces este archivo
             // mientras compila cierto .cpp

#define _USE_MATH_DEFINES

#include "mathematics.h"

///Some universal constants in cgs units

///Fundamental constants

const double cLight    = 2.9979e10;                   //speed of light (cm/s)
const double cLight2   = P2(cLight); 

const double planck    = 6.6261e-27;                  //Planck's constant (erg s)

const double pi        = 3.141592;

const double thomson   = 6.65e-25;                    //Thomson cross section

const double boltzmann = 1.3807e-16;			      //Boltzmann constant (erg/K)

const double fineStructConst		=  1.0/137.0;     //Fine structure constant

const double gravitationalConstant  =  6.674e-8;      //Gravitational constant (cm3 s-2 g-1)

const double stefanBoltzmann        = 5.6704e-5;      //Stefan-Boltzmann constant (erg/s*cm^2*K^4)	

const double bohrMagneton           = 9.27e-21;       //bohr magneton units?


///Astronomical Data

const double solarMass				=  1.998e33;               //Solar mass (g)
const double solarLuminosity		=  3.84e33;               //erg/s

///Unit conversion factors and others

const double EV_TO_ERG				= 1.602e-12;			   // eV to erg conversion

const double pc						= 3.09e18;			 	// parsec to cm

const double yr = 3.15e7;			 	// s in a year


///Particle Data 

const double electronCharge			=  4.8032e-10;             //electron charge (esu)

const double electronRadius			=  2.8179e-13;             //Electron classical                                                               //radius (cm)

//Particle rest mass energy (erg) mc^2
																
const double electronMass			= 9.1094e-28;   //(0.511e6);    
const double muonMass				= 1.8817e-25;   //(105.7e6);  
const double protonMass				= 1.6726e-24;   //(938.3e6);
const double chargedPionMass		= 2.4852e-25;   //(139.6e6);
const double neutralPionMass		= 2.4034e-25;   //(135.0e6);
const double neutronMass	     	= 1.6749e-24;   //(939.6e6);

//Particle mean life  (s)
                
const double chargedPionMeanLife	=  2.6e-8;	//s			   
const double muonMeanLife			=  2.2e-6;  //s	
const double neutronMeanLife        =  881.5;  //s

//Threshold energy

const double pionThresholdPH	=  200.0e6*1.6e-12;	 //145  pongo 200 y no 145 porque a partir de esta 
                                                     // energia esta parametrizada la seccion eficaz
const double pairThresholdPH	=  2.0*electronMass*cLight2;
const double pionThresholdH	    =  1.22e9*1.6e-12;
