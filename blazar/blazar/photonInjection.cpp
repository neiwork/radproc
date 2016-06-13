#include "photonInjection.h"

#include <fparameters\parameters.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>

#include <boost/property_tree/ptree.hpp>

#include <fluminosities\luminositySynchrotron.h>
#include <fluminosities\luminosityIC.h>
#include <fparameters/parameters.h>


double tpf(double E, const ParamSpaceValues psv, const SpaceCoord& distCoord)
{
	double result = psv.interpolate({ { 0, E } }, &distCoord);
	return result;
}



void photonTarget(Particle& p, State& st)
{
//	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double zInt = GlobalConfig.get<double>("zInt");

	double B = st.magf.get({ 0 });
	double Emin = p.emin();
	double Emax = eEmax(zInt, B);

	//volumen 
	double vol = pi*P2(jetRadius(zInt, openingAngle))*zInt;


	p.injection.fill([&](const SpaceIterator& i){

		double E = i.val(DIM_E);

		double Qsyn = luminositySynchrotron(E, st.electron, i, st.magf);

		double total = Qsyn*vol;

		return total;
	});

}

void photonDistribution(Particle& p, State& st)
{
	//	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double zInt = GlobalConfig.get<double>("zInt");

	double B = st.magf.get({ 0 });
	double Emin = p.emin();
	double Emax = eEmax(zInt, B);

	//volumen 
	double vol = pi*P2(jetRadius(zInt, openingAngle))*zInt;


	p.injection.fill([&](const SpaceIterator& i){

		double E = i.val(DIM_E);

		double Qsyn = luminositySynchrotron(E, st.electron, i, st.magf);

		const ParamSpaceValues psv = st.photon.injection;

		double Qic = luminosityIC(E, st.electron, i,
			[&psv, &i](double E){return tpf(E, psv, i); }
		, st.photon.emin());

		double total = Qsyn*vol;

		return total;
	});

}