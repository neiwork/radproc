#include "distribution.h"


//#include "injection.h"
#include "losses.h"

#include "modelParameters.h"
#include "messages.h"

#include <fmath\physics.h>
#include <flosses\lossesSyn.h>
#include <flosses\nonThermalLosses.h>
#include <fmath\RungeKutta.h>
#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\parameters.h>
#include <boost\property_tree\ptree.hpp>

void distribution(Particle& p, State& st)
{
	show_message(msgStart, Module_Distribution);

	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double zInt = GlobalConfig.get<double>("Rdiss");
	static const double gamma = GlobalConfig.get<double>("Gamma");

	const ParamSpace& ps{ p.ps };

	double v_lat = cLight*openingAngle;
	double B = st.magf.get({ 0 });
	double Emax = eEmax(p.mass, zInt, B);
	

		p.ps.iterate([&](const SpaceIterator& i){
			const double E = i.val(DIM_E);
				
			
			//double inj = p.injection.get(i);
			//double tad = E / adiabaticLosses(E, zInt, v_lat, gamma); //en [seg]

			//double tad =   (3.0*zInt)/ (2.0*cLight);

			//double tsyn = E / lossesSyn(E, B, p);

			//double tcross = zInt / cLight;

			//double tcool = 1.0/(1.0/tad + 1.0/tsyn);

			double integral = RungeKuttaSimple(E, Emax, [&](double e) {
				return p.injection.interpolate({ { DIM_E, E } });
			});

			double Ne = integral / (losses(E, p, st, i)); //Syn(E, B, p) + adiabaticLosses(E, zInt, v_lat, gamma));
			//double Qt = inj*tcool;
			p.distribution.set(i,Ne); 
				
		} );


}





