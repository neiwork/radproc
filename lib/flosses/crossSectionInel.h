#pragma once

/* Cross section and inelasticity for photopair interaction in the channel of pion production*/ 
double crossSectionPHPion(double t);
double inelasticityPHPion(double t);


/* Bethe Heitler cross section and inelasticity for photopair interaction in the channel of pair creation
	@param x is the normalized energy to the rest energy of electron*/ 
double crossSectionBetheHeitler(double t);
double inelasticityBetheHeitler(double t);


/*Cross section for pion production in pp/pip interactions*/ 
double crossSectionHadronic(double E);		// Kafexhiu et al. (2014)

/*Cross section for pion production in pp/pip interactions in delta aproximation*/ 
double crossSectionHadronicDelta(double Ekin);


double crossSectionThomson(double mass);

/*Klein-Nishina cross section*/ 
double crossSectionKN(double Eph, double Ee);
double angleAveragedKN(double eps); //eps is the normalized photon energy


/*Cross section for photon-photon pair production for frontal collision (theta = pi)
Tesis Mar'ia, Pe'er & Waxman 2004  */ 
double crossSectionGammaGamma(double beta);