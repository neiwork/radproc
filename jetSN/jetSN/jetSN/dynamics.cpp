#include "targetFields.h"

#include "modelParameters.h"

#include <iostream>
//#include <fmath\RungeKutta.h>
#include "State.h"
#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>
#include <fmath\physics.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>

double Fe(double g, double y)
{
	return pow(g, 4.0) * (1.0 / P2(g) - P2(g)) / P2(y);
}

//void gammaC(double z, double Gc)
void gammaC(State& st, Vector& Gc)
{
	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double theta = GlobalConfig.get<double>("openingAngle");
	static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double E_0 = GlobalConfig.get<double>("E_0");


	double Mc = 2.0*E_0 / P2(cLight);

	double z_int = st.electron.ps[DIM_R].first();

	double rc = stagnationPoint(z_int);

	double D = Lj*P2(rc) / (4.0*P2(theta)*P3(Gj*cLight)*z_int*Mc); //ver 

	std::cout << D;
	
	st.electron.ps.iterate([&](const SpaceIterator& i) {
		
		const double z = i.val(DIM_R);

		double y = z / z_int;
		
		int s = i.coord[DIM_R]; //posicion en la coordenada z
		double g;

		if (s == 0) {
			Gc[0] = 1.0;
			g = Gc[0] / Gj;
		}
		else {
			g = Gc[s - 1] / Gj; // 1.0e-2; comienzo la iteracion con Gc el que entra

			double y2 = 1.0;
			double err = 1.0;
			
			if (g < 1.0) {
				while (err > 1.0e-3 && y2 > 0.0) {

					y2 = pow((1.0 - (log((g + 1.0) / std::abs(g - 1.0)) - 2.0*std::atan(g)) / (4.0*D)), -1.0);

					g = g + 1.0e-3;
					Gc[s] = g*Gj;

					err = abs(y - y2) / y;

				}
			}
			else {
				Gc[s] = Gj;
			}
		}
	}, { 0, -1}); //fijo cualquier energia
}




double expansionTime(double z, double z_int)
{
	static const double E_0 = GlobalConfig.get<double>("E_0");

	double rc = stagnationPoint(z_int);
	double Mc = 2.0*E_0 / P2(cLight);

	double adiabatic = 4.0 / 3.0;

	double cs = sqrt(adiabatic*jetRamPress(z)*4.0*pi*rc / (3.0*Mc));

	double expT = 5.0*2.0*rc / cs; //A=5.0

	return expT;
}