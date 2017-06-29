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


void blobRadius(State& st, Vector& Gc, Vector& Rc)
{
	int nR = st.electron.ps[DIM_R].size();

	const double z0 = st.electron.ps[DIM_R].first();

	double R0 = stagnationPoint(z0);
	double cs = soundC(R0, z0);
	double t0 = z0 / cLight;
	
	st.electron.ps.iterate([&](const SpaceIterator& i) {

		const double z = i.val(DIM_R);
		int z_ix = i.coord[DIM_R]; //posicion en la coordenada z

		double gamma = Gc[z_ix];

		if (z_ix == 0) {
			Rc[z_ix] = R0;
		}
		else {
			/*Vector Raux(z_ix+1, 0.0);
			Raux[0] = Rc[0];
			for (size_t j = 1; j <= z_ix; ++j)
			{	
				double z_p = st.electron.ps[DIM_R][j - 1];
				double R_p = Raux[j - 1];
				double dt = (z - z_p) / cLight;

				Raux[j] = R_p + soundC(R_p, z_p)*dt / 10.0; // gamma;
			}*/
			
			double t = z / cLight - t0;

			double beta = cs / cLight;

			double R = R0 + cs*t / Gc[z_ix];
			
			Rc[z_ix] = R;
		}
	}, { 0, -1 }); //fijo cualquier energia

}



void gammaC(State& st, Vector& Gc, Vector& Rc)
{
	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double theta = GlobalConfig.get<double>("openingAngle");
	static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double E_0 = GlobalConfig.get<double>("E_0");


	double Mc = 2.0*E_0 / P2(cLight);

	double z0 = st.electron.ps[DIM_R].first();

	double R0 = stagnationPoint(z0);//aca dejo el R0, no el que se expande

	double cs = soundC(R0, z0);
	double t0 = z0 / cLight;


	//double cs = soundC(rc, z_int);
	//double texp = 5.0*R0 / cs;
	//double dt = z_int / (2.0*P2(Gj)*cLight*D);
	//double rmax = dt*cLight / pc;


	//std::cout << D << '\t' << rmax << '\t' << rc/jetRadius(z_int,theta) ;

	st.electron.ps.iterate([&](const SpaceIterator& i) {

		const double z = i.val(DIM_R);

		double y = z / z0;

		int s = i.coord[DIM_R]; //posicion en la coordenada z
		double g;


		if (s == 0) {
			Rc[0] = R0;
			Gc[0] = 1.0;
			g = Gc[0] / Gj;
		}
		else {

			//primero calculo Gc[s]
			double rc = Rc[s - 1];
			Rc[s] = rc;


			double D = Lj*P2(rc) / (4.0*P2(theta)*P3(Gj*cLight)*z0*Mc); //ver 

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

			//despues calculo Rc[s]

			double Plat = 1.0e-3*jetRamPress(z);
			double rhoj = Lj / (P2(Gj)*P3(cLight)*pi*P2(jetRadius(z, theta)));
			//P2(Gj)
			double vc = sqrt(1.0 - 1.0 / P2(Gc[s]));
			double vjet = sqrt(1.0 - 1.0 / P2(Gj));
			double v_rel = cLight*(vjet - vc);
			double G_rel = 1.0 / sqrt(1.0 - P2(v_rel / cLight));
			double Psn = rhoj*P2(v_rel*Gc[s]);
				
				//E_0 / (4.0*pi*P3(Rc[s]));

			double t = z / cLight - t0;

			double D2 = t*v_rel / z;

			if (Plat < Psn) {				

				Rc[s] = R0 + cs*t / Gc[s];
			}

			

		}
	}, { 0, -1 }); //fijo cualquier energia
}




double Fe(double g, double y)
{
	return pow(g, 4.0) * (1.0 / P2(g) - P2(g)) / P2(y);
}

//void gammaC(double z, double Gc)
/*void gammaC(State& st, Vector& Gc)
{
	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double theta = GlobalConfig.get<double>("openingAngle");
	static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double E_0 = GlobalConfig.get<double>("E_0");


	double Mc = 2.0*E_0 / P2(cLight);

	double z_int = st.electron.ps[DIM_R].first();

	double rc = stagnationPoint(z_int);//aca dejo el R0, no el que se expande

	double D = Lj*P2(rc) / (4.0*P2(theta)*P3(Gj*cLight)*z_int*Mc); //ver 


	//double cs = soundC(rc, z_int);
	//double texp = 5.0*R0 / cs;
	double dt = z_int / (2.0*P2(Gj)*cLight*D);
	double rmax = dt*cLight / pc;


	std::cout << D << '\t' << rmax << '\t' << rc/jetRadius(z_int,theta) ;
	
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
*/

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