#include "thermalSync.h"


#include <fmath/mathFunctions.h>
#include <fmath/physics.h>
#include <boost/math/special_functions/bessel.hpp>
#include "blackBody.h"
#include <gsl/gsl_sf_bessel.h>

double auxiliaryFunction(double& alpha, double& beta, double& gamma, double temp) {
    
    // se podría interpolar entre uno y otro.
    if (temp < 1.0e9) {
        alpha = 1.121;
        beta   = -10.65;
        gamma = 9.169;
    } else if (temp >= 1.0e9 && temp < 2.0e9) {
        alpha = 1.180;
        beta   = -4.008;
        gamma = 1.559;
    } else if (temp >= 2.0e9 && temp < 4.0e9) {
        alpha = 1.045;
        beta   = -0.1897;
        gamma = 0.0595;
    } else if (temp >= 4.0e9 && temp < 8.0e9) {
        alpha = 0.9774;
        beta   = 1.16;
        gamma = 0.2641;
    } else if (temp >= 8.0e9 && temp < 1.6e10) {
        alpha = 0.9768;
        beta   = 1.095;
        gamma = 0.8332;
    } else if (temp >= 1.6e10 && temp < 3.2e10) {
        alpha = 0.9788;
        beta   = 1.021;
        gamma = 1.031;
    } else if (temp >= 3.2e10) {
        alpha = 1.0;
        beta   = 1.0;
        gamma = 1.0;
    }
}

double mAux(double xM, double temp) {
    
    double alpha, beta, gamma;
    auxiliaryFunction(alpha, beta, gamma, temp);
	double result = (4.0505 * alpha / pow(xM, 1.0/6.0) ) * (1.0 + 0.4*beta / pow(xM, 0.25) + 
                0.5316 * gamma / sqrt(xM) ) * exp(-1.8899 * pow(xM, 1.0/3.0));
	    
    return result;
}

double jSync(double energy, double temp, double magfield, double dens_e)
{
	double frequency = energy / planck;
    double norm_temp = boltzmann * temp / electronRestEnergy;
	
    double nu0 = (electronCharge * magfield) / (2.0 *pi * electronMass * cLight);
    double xM = (2.0 * frequency) / (3.0 * nu0 * norm_temp*norm_temp);
	
	//double bessel = boost::math::cyl_bessel_k(2, 1.0/norm_temp);
	double bessel = gsl_sf_bessel_Kn(2,1.0/norm_temp);
	
    double result = (bessel > 0.0) ? 4.437e-30 * dens_e * frequency / bessel * mAux(xM, temp) : 0.0;
	return result;
} // esto debería tener unidades de erg cm^-3 ster^-1 s^-1 Hz^-1


/*double jSync(double energy, double temp, double magfield, double dens_e)
{
	double frecuency = energy / planck;
    double norm_temp = boltzmann * temp / (electronMass * cLight2);
	
    double nu0 = (electronCharge * magfield) / (2.0 *pi * electronMass * cLight);
    double xM = (2.0 * frecuency) / (3.0 * nu0 * norm_temp*norm_temp);
	
	//double bessel = boost::math::cyl_bessel_k(2, 1.0/norm_temp);
	double bessel = gsl_sf_bessel_Kn(2,1.0/norm_temp);
	
    double result = (bessel > 0.0) ? electronCharge*electronCharge / (sqrt(3.0)*cLight) *
	dens_e * frecuency / bessel * mAux(xM, temp) : 0.0;
    
	return result;
}*/ // esto debería tener unidades de erg cm^-3 ster^-1 s^-1 Hz^-1