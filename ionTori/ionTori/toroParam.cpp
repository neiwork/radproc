#include "toroParam.h"



#include <fparameters/SpaceIterator.h>
#include <fmath/bisection.h>
#include <fmath/physics.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>




// Keplerian specific angular momentum
double keplAngularMom(double r) {
		
	static const double massBH = GlobalConfig.get<double>("massBH");
    static const double spinBH = GlobalConfig.get<double>("spinBH") * massBH;

    return sqrt(massBH) * ( r*r - 2.0 * spinBH * sqrt(massBH*r) + spinBH*spinBH ) /
               (pow(r, 1.5) - 2.0 * massBH * sqrt(r) + spinBH* sqrt(massBH) );
}

// Torus Parameters
void torusParameters(double *l_0, double *rCusp, double *rCenter) {
    
	
	static const double massBH = GlobalConfig.get<double>("massBH");
    static const double spinBH = GlobalConfig.get<double>("spinBH") * massBH;
	static const double lambda = GlobalConfig.get<double>("lambda");
	
	// Auxiliary variables
	
    double z1 = 1.0 + pow( 1.0 - (spinBH/massBH)*(spinBH/massBH), 1.0/3.0) * 
    ( pow(1.0 + spinBH/massBH, 1.0/3.0) + pow(1.0 - spinBH/massBH, 1.0/3.0) );
    double z2 = sqrt( 3.0 * (spinBH/massBH)*(spinBH/massBH) + z1*z1);

    // marginally stable circular orbit
    double r_ms = massBH * (3.0 + z2 - sqrt( (3.0 - z1) * (3.0 + z1 + 2.0*z2) ) );

    // marginally bound circular orbit
    double r_mb = 2.0 * massBH - spinBH + 2.0 * sqrt(massBH) * sqrt(massBH-spinBH);

    double l_ms = keplAngularMom(r_ms);              // Keplerian specific angular momentum at r = r_ms
    double l_mb = keplAngularMom(r_mb);             // Keplerian specific angular momentum at r = r_mb
    
    *l_0 = (1.0 - lambda) * l_ms + lambda * l_mb;

	*rCusp = bisection(r_mb, r_ms, [&l_0](double r) {return keplAngularMom(r) - (*l_0); } );
	*rCenter = bisection(r_ms, 1000.0, [&l_0](double r) {return keplAngularMom(r) - (*l_0); } );

}
