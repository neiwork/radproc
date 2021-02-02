#include "ebl_absorption.h"

#include <iostream>
#include <fstream>

void absorption(Vector& E_TeV, Vector& tau02, Vector& tau1)
{
	std::ifstream fileReaded;
	fileReaded.open("ebl_tau.txt", std::ios::in);


	//double E_TeV, tau02, tau1;
	int j = 0;


	for(size_t j=0; j < E_TeV.size(); ++j)
	{
		fileReaded >> E_TeV[j] >> tau02[j] >> tau1[j];
	}

	fileReaded.close();

}
