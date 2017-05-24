#include "read.h"

#include <fparameters\parameters.h>

//using namespace std;

namespace {
	double safeLog10( double x ) {
		return x>0.0 ? log10(x) : 0.0;
	}
}

//
//void read (const std::string& archive, Vector& ener, Vector& lum)
//{
//	std::ifstream fileReaded;
//	fileReaded.open(archive.c_str(), std::ios::in);
//
//	double logE, loglum;
//	int j=0;
//
//	while (! fileReaded.eof())  //esto termina cuando llega al final
//		{			
//			fileReaded >> logE;
//			fileReaded >> loglum;
//			ener[j]  = pow(10.0,logE)*1.6e-12;
//			lum[j] = pow(10.0,loglum);   
//			j += 1;
//		}
//
//	fileReaded.close();	
//}
//
//void readChanged (const std::string& archive, Vector& ener, Vector& lum)
//{
//	std::ifstream fileReaded;
//	fileReaded.open(archive.c_str(), std::ios::in);
//
//	double logE, loglum;
//	int j=0;
//
//	while (! fileReaded.eof())  //esto termina cuando llega al final
//		{			
//			fileReaded >> logE;
//			fileReaded >> loglum;
//			ener[j]  = pow(10.0,logE)*1.6e-12;
//			lum[j] = pow(10.0,loglum)*radius/(volume*P2(ener[j])*cLight);  
//			j += 1;
//		}
//
//	fileReaded.close();	
//}

