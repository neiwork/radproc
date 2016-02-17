#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <fparticle\Particle.h>
#include <fmath\physics.h>

class File;

typedef std::map<std::string, File*> OFM;

double safeLog10(double x);

class File {
public:
	std::string name;
	std::ofstream file;
	File(std::vector<File*>& all, OFM& byName, std::string name) :name(name) {
		file.open(name + ".txt");
		byName[name] = this;
		all.push_back(this);
	}
	void check() {
		std::ifstream result;
		result.open(name + ".txt");

		std::ifstream reference;
		reference.open(name + ".ref");

		char a, b;
		int c = 0;
		while (!result.eof()) {
			c++;
			result >> a;
			reference >> b;
			if (a != b) {
				std::cout << name << +" doesn't match reference." << std::endl;
				throw;
			}
		}
		if (!reference.eof()) {
			std::cout << name << +" reference content larger than result." << std::endl;
			throw;
		}
		result.close();
		reference.close();
		std::cout << name << +" matches reference (" << c << "b)." << std::endl;
	}
};

void writeFlux(const std::string& archive, Vector salida, Particle particle);

void write(const std::string& archive, Vector salida, Particle particle);

void writeAllSpaceParam(const std::string& filename, const ParamSpaceValues& data);

void writeEnergyFunction(const std::string& filename, const ParamSpaceValues& data, int r, int t);

/* Write is a function that given the vectors ener and lum, prints on archive logE and logL.  */
//void write(const std::string& archive, Vector& ener, Vector& salida);

/* WriteChanged is a function that given the vectors ener and lum,
prints on archive logE and logL with some changes in the midddle.  */
//void writeChanged(const std::string& archive, Vector& ener, Vector& lum);

/*Function writeMatrix writes a particular matrix*/
//void writeMatrix (const std::string& archive, Matrix external, Particle particle);

/*Function writeLum writes the total luminosity*/
//void writeLum(const std::string& archive, Vector salida, Particle particle);



//void writeEvolution (const std::string& electronFilename, 
//					 const std::string& positronFilename, 
//					 const std::string& photonFilename, CESEvolution states);

//void writeMatrix(const std::string& archive, Matrix external, Particle particle);