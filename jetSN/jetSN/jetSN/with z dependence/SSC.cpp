#include "SSC.h"

#include "processes.h"

#include <fluminosities\luminositySynchrotron.h>
#include <fluminosities\luminosityIC.h>

#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>

#include <fparameters/parameters.h>

#include <fmath\physics.h>

#include <boost/property_tree/ptree.hpp>




double nSSC(double E, ParamSpace pps, ParamSpaceValues Qsyn)
{
	int t_ix = pps[DIM_R].size() - 1;

	double t = pps[2][t_ix];  //t es el ultimo valor

	double Lsyn = emiToLumi(pps, Qsyn, E, t_ix);
	return Lsyn;

}


void SSC(State st, ParamSpaceValues& Qssc, ParamSpaceValues Qsyn)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	double r = st.electron.ps[DIM_R].first();

	const ParamSpace& pps = st.photon.ps;

	st.photon.ps.iterate([&](const SpaceIterator &i){

		const double E = i.val(DIM_E);		

		const double eSSC = luminosityIC(E, st.electron, i,
			[&Qsyn, &pps, &r](double E){return nSSC(E, pps, Qsyn) / (P2(E) *4.0*pi*P2(jetRadius(r, openingAngle))*cLight); }
		, st.photon.emin());

		//const double eSSC = luminosityIC(E, st.electron, i.coord, [&](double E){
		//	return luminositySynchrotron(E, st.electron, i.coord, st.magf) /
		//		(P2(E) *4.0*pi*P2(jetRadius(r, openingAngle))*cLight); }, st.photon.emin());
		Qssc.set(i, eSSC);
	}, { -1, -1, (int)st.photon.ps[DIM_R].size() - 1 });
}







//double tpf(double E, const ParamSpaceValues psv, const SpaceCoord& distCoord)
//{
//	double result = psv.interpolate({ { 0, E } }, &distCoord);
//	return result;
//}