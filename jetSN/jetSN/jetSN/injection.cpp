#include "injection.h"

#include "messages.h"
#include "modelParameters.h"
#include "nonThermalLuminosity.h"

#include <fparameters\parameters.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\Dimension.h>

#include <fmath\RungeKutta.h>
#include <fmath\physics.h>

#include <boost/property_tree/ptree.hpp>

#include <iostream>

double powerLaw(double E, double Emin, double Emax)
{
	static const double primaryIndex = GlobalConfig.get<double>("primaryIndex");
	
	double result = pow(E, (-primaryIndex))*exp(-E / Emax)*exp(-5 * Emin / E);
	return result;
}


double normalization(Particle& p, double z, double B)
{
	static const double Gamma = GlobalConfig.get<double>("Gamma");

	//int i_z = distCoord[1];
	//double z = p.ps[1][i_z];

	double Emin = p.emin();
	double Emax = eEmax(z, B);

	double int_E = RungeKuttaSimple(Emin, Emax, [&Emax, &Emin](double E){
		return E*powerLaw(E, Emin, Emax);
	});  //integra E*Q(E)  entre Emin y Emax
	
	double Q0 = dLnt(z) / (int_E);  //factor de normalizacion de la inyeccion
	return Q0 / P2(Gamma);
	//N'(E')dE' = N(E)dE  ==> E'_nt = int( E'N(E')dE') = E_nt/Gamma;
}



void injection(Particle& p, State& st)
{
	//static const double magf = GlobalConfig.get<double>("magf");

	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double z_int = GlobalConfig.get<double>("z_int")*pc;

	show_message(msgStart, Module_electronInjection);

	//const double RMIN = p.ps[DIM_R].first();
	//const double RMAX = p.ps[DIM_R].last();
	//const int N_R = p.ps[DIM_R].size()-1;

	double Lnt_total = nonThermalLuminosity();

	//volumen total del jet
	//double vol = (pi / 3.0)*(P2(jetRadius(RMAX, openingAngle))*RMAX
	//	- P2(jetRadius(RMIN, openingAngle))*RMIN);

	//double z_int = pow((RMAX / RMIN), (1.0 / N_R));

	//int z_ix = 0;
	//int t_ix = 0;

	//static const std::string injector = GlobalConfig.get<std::string>("injector");

	//bool multiple = (injector == "multiple");
	//bool single = (injector == "single");
	//bool condicion;

	

	p.injection.fill([&](const SpaceIterator& i){
		const double magf{ st.magf.get(i) };
		//const double r{ i.val(DIM_R) };
		
		//if (single)
		//{
		//	condicion = i.its[2].canPeek(-1) || i.its[1].canPeek(-1);
			//if (i.its[2].canPeek(-1) || i.its[1].canPeek(-1)) /* injector en z=0 */	
		//}
		//else if (multiple)
		//{
		//	condicion = i.its[2].canPeek(-1);
		//	//if (i.its[2].canPeek(-1))     /* injectores para todo z */
		//}

		//if (condicion) /* injector en z=0 */
		//{
		//	return 0.0;
		//}
		//else //if (t_position = 0) solo inyecto particulas a tiempo 0
		//{
		
			//if (i.its[1].dim == 0 && i.its[2].dim){
			double Emin = p.emin();
			double Emax = eEmax(z_int,magf);
			double Q0 = normalization(p,z_int,magf);
		
			//double z = i.val(DIM_R);			
			//double dz = z*(z_int - 1);
			//volumen de la celda i
			//double vol_i = pi*P2(jetRadius(z, openingAngle))*dz;

			double total = powerLaw(i.val(DIM_E), Emin, Emax)*Q0; // *vol_i;// / vol;

			return total;
		//}

	});

	show_message(msgEnd, Module_electronInjection);
}


