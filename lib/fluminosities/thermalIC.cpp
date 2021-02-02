#include "thermalIC.h"

#include <fmath/physics.h>
#include <fmath/RungeKutta.h>
#include <fparameters/parameters.h>
#include <fparameters/ParamSpace.h>
#include <fparameters/ParamSpaceValues.h>
#include <iostream>
#include <fmath/fbisection.h>

void auxFunction(double& e1, double& e2,  double& e3, double xc, double energy, double norm_temp, 
                                double r, double theta, const SpaceCoord& distCoord, ParamSpaceValues& denf) {
    
	double rg = gravitationalConstant * 1.0e6 * solarMass / cLight2;
    
	r = r / rg;
	double lim_inf = r*1.01;
    double lim_sup = r*10.0;
 
    double opticalDepth = rg*RungeKuttaSimple(lim_inf, lim_sup, [&](double r){
		try {
			return thomson * denf.interpolate({ {DIM_R, r} , {DIM_THETA, theta} }, &distCoord);
		} catch (std::runtime_error& e) {
			std::cout << "WARNING: " << e.what() << std::endl;
			return 0.0;
		}
	});
                    
    double P = 1.0 - exp(-opticalDepth);
    double A = 1.0 + norm_temp * (4.0 + 16.0 * norm_temp);
    e1 = P * (A - 1.0) / (1.0 - A*P);
    e3 = -1.0 - log(P) / log(A);
    e2 = pow(3.0, -e3) * e1;
}

double jIC_Bremss(double energy, double norm_temp, double r,  double theta, const SpaceCoord& distCoord, 
						ParamSpaceValues& denf, double jBr, double xc) {

    double eta1 = 0.0, eta2 = 0.0, eta3 = 0.0;
	double jBremss = jBr*energy*energy * 0.25 / pi; //porque el tpf esta /E^2
    
    double x = energy / (electronMass*cLight2);
	
	if (x < 3.0 * norm_temp)
	{
		auxFunction(eta1, eta2, eta3, x, energy, norm_temp, r, theta, distCoord, denf);
		
		double coeff1 = jBremss*eta1*norm_temp;
		double coeff2 = ( 1.0 - xc/norm_temp);
		double coeff3 = (3.0/(eta3+1.0))* ( pow(3.0,-eta3-1.0) - pow(x/(3.0*norm_temp) , eta3+1.0) ) ;
		double result = coeff1* ( coeff2 - coeff3);
	
		return result;
		
	}
	else {	
		return 0.0;
		} //> 0.0 ? result : 0.0;
}
    
double jIC_Sync(double energy, double norm_temp, double r, double theta, const SpaceCoord& distCoord,
                ParamSpaceValues& denf, double jSync, double xc) {
     
    double x = energy / (electronMass * cLight2);
    
    if (x < 3.0 * norm_temp) {
        double eta1 = 0.0, eta2 = 0.0, eta3 = 0.0, xc = 0.0;
        auxFunction(eta1, eta2, eta3, x, energy, norm_temp, r, theta,
                            distCoord, denf);
    
        return jSync * ( eta1 - eta2*pow(x/norm_temp, eta3) );
    } else {
        return 0.0;
    }
}

double jIC(double jSource, double normTemp, double r, double rMax, double theta, const SpaceCoord& distCoord,
                ParamSpaceValues& denf, int s) {
    
    double rg = gravitationalConstant * 1.0e6 * solarMass / cLight2;
    
	//r = r / rg;
	double lim_inf = r;//*1.01;
    double lim_sup = rMax;//*0.9;
 
    double opticalDepth = rg*RungeKuttaSimple(lim_inf, lim_sup, [&](double r){
		try {
			return thomson * denf.interpolate({ {DIM_R, r} , {DIM_THETA, theta} }, &distCoord);
		} catch (std::runtime_error& e) {
			std::cout << "WARNING: " << e.what() << std::endl;
			std::cout << lim_inf << '\t' << lim_sup << '\t' << r << std::endl;
			return 0.0;
		}
	});
    
    double probEscape = exp(-opticalDepth);
    double probScattering = 1.0 - probEscape;
    
    double probability = probEscape;
    for (int i=0; i < s; i++) {
        probability *= probScattering;
    }
    
    return jSource * probability;
}