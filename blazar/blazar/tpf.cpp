#include "tpf.h"

#include "modelParameters.h"

#include <fparameters\SpaceIterator.h>
#include <fparameters\parameters.h>
#include <fluminosities\luminositySynchrotron.h>

//#include <boost/property_tree/ptree.hpp>
//#include <fparameters\Dimension.h>



/* Takes [emi] =  E^2*[Q(E)] and calculates int(2.0*pi*P2(jetR)*emi dz);
for [N(E)] = 1/erg, then it just sums over all z and returns erg/s  */

void synL(State& st, ParamSpaceValues& Qsyn)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double zInt = GlobalConfig.get<double>("Rdiss");

	st.photon.ps.iterate([&](const SpaceIterator &i) {

		const double E = i.val(0);
		double eSyn = luminositySynchrotron(E, st.electron, i.coord, st.magf); //estos devuelven erg/s, sumar!

		Qsyn.set(i, eSyn / (P2(E) *4.0*pi*P2(jetRadius(zInt, openingAngle))*cLight));
	});
}




void photonTarget(Particle& p, State& st)
{
	//	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double zInt = GlobalConfig.get<double>("Rdiss");


	double magf = st.magf.get({ 0 });
	//double Emin = p.emin();
	//double Emax = p.emax();//es la del foton eEmax(zInt, B);

	double radius = jetRadius(zInt, openingAngle);
	//double vol = pi*P2(radius)*zInt;

	p.injection.fill([&](const SpaceIterator& i) {

		double E = i.val(DIM_E);

		double Qsyn = luminositySynchrotron(E, st.electron, i, st.magf); // st.magf); //erg/s  double luminositySynchrotron(double E, const Particle& c, const SpaceCoord& psc, double magf)

		double total = Qsyn / (P2(E) *4.0*pi*P2(radius)*cLight); //lo paso a 1/erg cm^3

		return total;
	});

}