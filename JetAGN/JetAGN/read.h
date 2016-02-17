#pragma once

//#include "CoupledEqSys.h"
#include <string>
#include <fparticle\Particle.h>
#include <iostream>
#include <fstream>
#include <fmath\physics.h>




/* Read is a function that from archive obtained logE and logL and generates the vectors ener and lum. */ 
void read (const std::string& archive, Vector& ener, Vector& lum);

/* ReadChanged is a function that from archive
obtained logE and logL and generates the vectors ener and lum with some changes in the midddle. */ 
void readChanged (const std::string& archive, Vector& ener, Vector& lum);
