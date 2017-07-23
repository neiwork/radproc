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



double soundC(double R, double z)
{
	static const double E_0 = GlobalConfig.get<double>("E_0");
	const double adiaCoeff = 4.0 / 3.0;
	double Mc = 2.0*E_0 / cLight2;

	double cs = sqrt(4.0*pi*P3(R)*adiaCoeff*jetRamPress(z) / (3.0*Mc));

	double beta = cs / cLight;
	return cs;
}



void gammaC(State& st, Vector& Gc, Vector& Rc, Vector& tobs)
{
	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double theta = GlobalConfig.get<double>("openingAngle");
	static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double E_0 = GlobalConfig.get<double>("E_0");
	static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree
	
	double mu = cos(inc);

	const double z0 = st.electron.ps[DIM_R].first();
	const double zMax = st.electron.ps[DIM_R].last();
	const int N_R = st.electron.ps[DIM_R].size() - 1;

	double z_int = pow((zMax / z0), (1.0 / N_R));
	
	double t0 = z0 / cLight;
	double t = 0.0;
	
	double Mc = 2.0*E_0 / P2(cLight);
	
	double R0 = stagnationPoint(z0);//aca dejo el R0, no el que se expande

	double cs = soundC(R0, z0);
	//double cs = soundC(rc, z_int);

	st.electron.ps.iterate([&](const SpaceIterator& i) {

		const double z = i.val(DIM_R);

		double y = z / z0;

		int s = i.coord[DIM_R]; //posicion en la coordenada z
		double g;
		
		if (s == 0) {
			Rc[0] = R0;
			Gc[0] = 1.0;
			g = Gc[0] / Gj;
			tobs[0] = 0.0;
		}
		else {

			//primero calculo Gc[s]
			double rc = Rc[s - 1];
			Rc[s] = rc;
			
			double D = Lj*P2(rc) / (4.0*P2(theta)*P3(Gj*cLight)*z0*Mc); 

			g = Gc[s - 1] / Gj; // 1.0e-2; comienzo la iteracion con Gc el que entra

			double y2 = 1.0;
			double err = 1.0;

			if (g < 1.0) {
				while (err > 1.0e-3 && y2 > 0.0) {

					//y2 = pow((1.0 - (log((g + 1.0) / std::abs(g - 1.0)) - 
					//     2.0*std::atan(g)) / (4.0*D)), -1.0);  //mhd

					y2 = pow((1.0 - (log(std::abs(g - 1.0)/(g + 1.0)) - 
						2.0*g/(P2(g)-1.0)) / (4.0*D*Gj)), -1.0); //hydro

					g = g + 1.0e-3;
					Gc[s] = g*Gj;

					err = abs(y - y2) / y;

				}
			}
			else {
				Gc[s] = Gj;
			}

			//Calculo Rc[s]

			double beta_c = sqrt(1.0 - 1.0 / P2(Gc[s]));
			double beta_j = sqrt(1.0 - 1.0 / P2(Gj));

			double beta_rel = (beta_j - beta_c) / (1.0 - beta_j*beta_c);
			double G_rel = 1.0 / sqrt(1.0 - P2(beta_rel));

			double Psn = P2(beta_rel*G_rel);
			double Plat = 1.0e-3*P2(beta_j*Gj);

			double dt = z*(z_int - 1.0) / (beta_c*cLight);   // agregue el beta_c;

			if (Plat < Psn) {

				Rc[s] = Rc[s-1] + cs*dt / Gc[s];
			}

			double cte = (1.0 - beta_c*mu);
			tobs[s] = tobs[s - 1] + dt*cte;
		}

		//}
	}, { 0, -1 }); //fijo cualquier energia
}




double Fe(double g, double y)
{
	return pow(g, 4.0) * (1.0 / P2(g) - P2(g)) / P2(y);
}


double eEmax(double z0, double z, double Gc, double B, double Reff)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");

	//double Reff = blobRadius(p, Gc, z_ix);
	//double Reff = blobRadius(z0, z, Gc);
	double vel_lat = cLight*openingAngle;

	double Emax_ad = accEfficiency*3.0*jetRadius(z, openingAngle)*cLight*electronCharge*B / (vel_lat*Gc); //
	double Emax_syn = electronMass*cLight2*sqrt(accEfficiency*6.0*pi*electronCharge / (thomson*B));

	double ampl = Gc; //factor de amplificaci'on de B en la zona del choque
	double Emax_diff = electronCharge*B*Reff*sqrt(3.0*accEfficiency*ampl / 2.0);
	double min1 = std::min(Emax_syn, Emax_ad);
	double min2 = std::min(min1, Emax_diff);


	return min2;

}

void fillMagnetic(State& st, Vector& Gc)
{
	static const double Gj = GlobalConfig.get<double>("Gamma");

	st.magf.fill([&](const SpaceIterator& i) {
		double z = i.val(DIM_R);
		int z_ix = i.coord[DIM_R];

		double beta_c = sqrt(1.0 - 1.0 / P2(Gc[z_ix]));
		double beta_j = sqrt(1.0 - 1.0 / P2(Gj));
		double beta_rel = (beta_j - beta_c) / (1.0 - beta_j*beta_c);
		double G_rel = 1.0 / sqrt(1.0 - P2(beta_rel));

		return computeMagField(z, G_rel);  // 
	});
}







/* 

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


void blobRadius(State& st, Vector& Gc, Vector& Rc)
{
	static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree
	double mu = cos(inc);

	const double z0 = st.electron.ps[DIM_R].first();
	const double zMax = st.electron.ps[DIM_R].last();
	const int N_R = st.electron.ps[DIM_R].size() - 1;

	double z_int = pow((zMax / z0), (1.0 / N_R));


	double t0 = z0 / cLight;
	double t = 0.0;

	double R0 = stagnationPoint(z0);
	double cs = soundC(R0, z0);


	st.electron.ps.iterate([&](const SpaceIterator& i) {

		const double z = i.val(DIM_R);
		int z_ix = i.coord[DIM_R]; //posicion en la coordenada z

		double gamma = Gc[z_ix];

		if (z_ix == 0) {
			Rc[z_ix] = R0;
		}
		else {
			

			double t = z / cLight - t0;

			double beta = cs / cLight;

			double R = R0 + cs*t / Gc[z_ix];

			Rc[z_ix] = R;
		}
	}, { 0, -1 }); //fijo cualquier energia

}

*/ 