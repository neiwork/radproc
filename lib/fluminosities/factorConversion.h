#pragma once


/* this factor is to convert from luminosity [erg/s] to Nph [1/erg cm^3]  */
double factorLumToNph(double E);


/* this factor is to convert from luminosity [erg/s] to Qph [1/erg cm^3 s]  */
double factorLumToQph(double E);