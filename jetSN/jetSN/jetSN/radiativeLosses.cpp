#include "radiativeLosses.h"

#include "dynamics.h"
#include "targetFields.h"
#include "write.h"
#include "messages.h"

#include <flosses\lossesSyn.h>
#include <flosses\nonThermalLosses.h>
#include <flosses\lossesIC.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>

void radiativeLosses(State& st, const std::string& filename, Vector& Gc, Vector& Rc)
{
	show_message(msgStart, Module_radLosses);

	static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");
	static const double starT = GlobalConfig.get<double>("starT");
	static const double starTIR = GlobalConfig.get<double>("IRstarT");
	static const double inc = GlobalConfig.get<double>("inc")*pi / 180;

	double z_0 = st.electron.ps[DIM_R].first();
	double t0 = z_0 / cLight;
	double vjet = sqrt(1.0 - 1.0 / P2(Gj));

	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << "Log(E/eV)" 
		<< "\t" << "z [pc]"
		<< "\t" << "tobs [yr]"
		<< "\t" << "Synchr"
		<< "\t" << "IC"
		<< "\t" << "IC - IR"
		<< "\t" << "Esc"
		<< "\t" << "Acc"
		//<< "\t" << "Ad"
		<< std::endl;

	double Emin = boltzmann*starT / 100.0;
	double EphminAux = boltzmann*starTIR / 100.0;


	st.electron.ps.iterate([&](const SpaceIterator& i) {
		const double B = st.magf.get(i);
		double vel_lat = cLight*openingAngle;

		double E = i.val(DIM_E);
		double z = i.val(DIM_R);

		double gamma = Gc[i.coord[DIM_R]];

		double Reff = z / gamma; // Rc[i.coord[DIM_R]];

		double fmtE = log10(E / 1.6e-12);

		double eSyn = lossesSyn(i.val(DIM_E), B, st.electron) / i.val(DIM_E);

		//VER como le paso el vector Gc a las perdidas
		double eIC = lossesIC(i.val(DIM_E), st.electron,
			[&E, &z, &gamma](double E) {
			return starBlackBody(E, z, gamma); },
			Emin, 1.0e4*Emin) / i.val(DIM_E);

			double eIC_Aux = lossesIC(i.val(DIM_E), st.electron,
				[&E, &z, &gamma](double E) {
				return starIR(E, z, gamma); },
				EphminAux, 1.0e4*EphminAux) / i.val(DIM_E);


				double vc = sqrt(1.0 - 1.0 / P2(gamma));
				double v_rel = cLight*(vjet - vc);
				double eEsc = escapeRate(Reff, v_rel);

				double eAcc = accelerationRate(E, B, accEfficiency);
				//double eAdia = adiabaticLosses(E, z, vel_lat, gamma) / E;

				double tlab = (z / cLight - t0) / yr;
				double t = tobs(tlab, gamma, inc);

				file << fmtE << '\t' << z / pc
					<< '\t' << t
					<< "\t" << safeLog10(eSyn)
					<< "\t" << safeLog10(eIC)
					<< "\t" << safeLog10(eIC_Aux)
					<< "\t" << safeLog10(eEsc)
					<< "\t" << safeLog10(eAcc)
					//<< "\t" << safeLog10(eAdia)
					<< std::endl;


	});


	file.close();
}


