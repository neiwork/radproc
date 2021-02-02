#pragma once
#include <../ionTori/modelParameters.h>

/* Straub et al 2012; [jIC] = erg cm^-3 ster^-1 s^-1 Hz^-1 */

double jIC_Bremss(double energy, double norm_temp, double r, double theta, const SpaceCoord& distCoord, 
                        ParamSpaceValues& denf, double jBr, double xc);

double jIC_Sync(double energy, double norm_temp, double r, double theta, const SpaceCoord& distCoord, 
                        ParamSpaceValues& denf, double jSy, double xc);
                        
double jIC(double jSource, double normTemp, double r, double rMax, double theta, 
                    const SpaceCoord& distCoord, ParamSpaceValues& denf, int s);
                        