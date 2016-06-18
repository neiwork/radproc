#include "distribution.h"


#include "injection.h"
//#include "losses.h"

#include "modelParameters.h"

#include <fmath\physics.h>
#include <flosses\lossesSyn.h>
#include <flosses\nonThermalLosses.h>
#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\parameters.h>
#include <boost/property_tree/ptree.hpp>

void distribution(Particle& p, State& st)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double zInt = GlobalConfig.get<double>("Rdiss");

	const ParamSpace& ps{ p.ps };

	double v_lat = cLight*openingAngle;
	double B = st.magf.get({ 0 });
	double Emax = eEmax(zInt, B);
	

		p.ps.iterate([&](const SpaceIterator& i){
			const double E = i.val(DIM_E);
				
			
			double inj = p.injection.get(i);

			double tad = E / adiabaticLosses(E, zInt, v_lat); //en [seg]
			double tsyn = E / lossesSyn(E, B, p);

			double tcool = std::min(tad, tsyn);
			p.distribution.set(i,inj*tcool); 
				
		} );


}





