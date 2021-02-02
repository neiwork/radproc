#include "sscLosses.h"

#include "lossesAnisotropicIC.h"

//#include "dynamics.h"
#//include "targetFields.h"
#include "write.h"
#include "messages.h"

//#include <flosses\lossesPhotoHadronic.h>
//#include <flosses\lossesSyn.h>
//#include <flosses\nonThermalLosses.h>
#include <flosses\lossesIC.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>

void sscLosses(State& st, const std::string& filename, Vector& Gc, Vector& tobs, ParamSpaceValues Qsyn)
{
	show_message(msgStart, Module_radLosses);

	//static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	//static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	//static const double starT = GlobalConfig.get<double>("starT");
	//static const double starTIR = GlobalConfig.get<double>("IRstarT");
	static const double z_peak = GlobalConfig.get<double>("z_peak")*pc;

	double z_0 = st.electron.ps[DIM_R].first();
	double t0 = z_0 / cLight;
	
	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << "Log(E/eV)"
		<< "\t" << "z [pc]"
		<< "\t" << "tobs [yr]"
		//<< "\t" << "Synchr"
		<< "\t" << "IC-ssc"
		//<< "\t" << "IC - IR"
		//<< "\t" << "Esc"
		//<< "\t" << "Acc"
		//<< "\t" << "Ad"
		<< std::endl;

	
	st.electron.ps.iterate([&](const SpaceIterator& i) {
		

		double E = i.val(DIM_E);
		double z = z_peak;  // i.val(DIM_R);

		double gamma = Gc[i.coord[DIM_R]];


		double fmtE = log10(E / 1.6e-12);

		double eIC =	lossesIC_old(i.val(DIM_E), st.electron,
			[&Qsyn, &i, &z, &st](double E) {
			//return Qsyn.interpolate({ { DIM_E, E },{ DIM_R, r } }) / (P2(E) *4.0*pi*P2(jetRadius(r, openingAngle))*cLight); }
			if (E < st.photon.emax() && E > st.photon.emin()) {
				return Qsyn.interpolate({ { DIM_E, E },{ DIM_R, z } }) / (P2(E) *4.0*pi*P2(jetRadius(z, openingAngle))*cLight);
			}
			else { return 0.0; }; }, 
			st.photon.emin(),st.photon.emax()) / i.val(DIM_E);



				double tlab = (z / cLight - t0) / yr;
				double t = tobs[i.coord[DIM_R]];

				file << fmtE << '\t' << z / pc
					<< '\t' << t / yr
					<< "\t" << safeLog10(eIC)
					<< std::endl;


	}, { -1, 0 });  //solo para el primer z


	file.close();
}


