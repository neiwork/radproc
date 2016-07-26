#include "SSC.h"

#include "processes.h"
#include "targetFields.h"

#include <fluminosities\luminositySynchrotron.h>
#include <fluminosities\luminosityIC.h>

#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>

#include <fparameters/parameters.h>

#include <fmath\physics.h>

#include <boost/property_tree/ptree.hpp>




double nSSC(double E, double r, ParamSpace pps, ParamSpaceValues Qsyn)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	int t_ix = pps[DIM_R].size() - 1;

	double t = pps[2][t_ix];  //t es el ultimo valor

	double Lsyn = emiToLumi(pps, Qsyn, E, t_ix);
	return Lsyn/(P2(E) *4.0*pi*P2(jetRadius(r, openingAngle))*cLight);

}

void SSC(State st, ParamSpaceValues& Qssc, ParamSpaceValues Qsyn)
{
	const ParamSpace& pps = st.photon.ps;

	ParamSpaceValues SynL(st.photon.ps);
	tpfPSV(SynL, 
		[&pps, &Qsyn](double E, double z){return nSSC(E, z, pps, Qsyn); }
		, st.photon, 1.0);

	st.photon.ps.iterate([&](const SpaceIterator &i){

		const double E = i.val(DIM_E);
		const double r = i.val(DIM_R);

		const double eSSC = luminosityIC(E, st.electron, i, SynL, st.photon.emin());
			//[&Qsyn, &pps, &r](double E){return nSSC(E, pps, Qsyn) / ; }
		

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