#include "luminosityPhotoHadronic.h"

#include <fparameters\parameters.h>
#include <fmath\RungeKutta.h>
#include <fmath\interpolation.h>
#include <flosses\crossSectionInel.h>
#include <finjection\pgammaPionInj.h>
#include <fmath\physics.h>

//ouble fHadron(double x, const Particle& creator,
	//const ParamSpaceValues& denf, const SpaceCoord& psc) //funcion a integrar   x=Ecreator; L=L(Ega)
double luminosityPhotoHadronic(double E, const Particle& creator, fun1 tpf, const SpaceCoord& psc, double tpEmin, double tpEmax)
{
	double diezE = 10.0*E;
	
	double distCreator = creator.distribution.interpolate({ { 0, diezE } }, &psc);

	double t_1   = t_pion_PH(diezE, creator, tpf, tpEmin, tpEmax);     //esto no es lossesPH porque son perdidas solo del canal de produccion de piones
	double omega = omegaPH(diezE, creator, tpf, tpEmin, tpEmax);

	if (omega > 0.0)	{
		double averageInel = t_1/omega;

		double k1 = 0.2;

		double k2 = 0.6;

		double p1 = (k2-averageInel)/(k2-k1);
	
		double nNeutralPion = 1.0-0.5*p1;

		double luminosity = 20.0*nNeutralPion*omega*distCreator;

		return luminosity*P2(E);  // [erg s^-1 cm^-3 ]
	}
	else	{return 0.0;}

	
}
