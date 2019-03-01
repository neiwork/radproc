#include "targetFields.h"

#include "modelParameters.h"

#include <iostream>
//#include <fmath\RungeKutta.h>
#include "State.h"
#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>
#include <fmath\bisection.h>
#include <fmath\physics.h>
#include <fparameters/parameters.h>

#include <string>
#include <iostream>
#include <fstream>

#include <boost/property_tree/ptree.hpp>



double soundC(double R, double z)
{
	static const double Lj = GlobalConfig.get<double>("Lj");
	static const double theta = GlobalConfig.get<double>("openingAngle");
	static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double g = GlobalConfig.get<double>("g");
	static const double Mc = GlobalConfig.get<double>("Mc")*solarMass;
		

	const double adiaCoeff = 4.0 / 3.0;
	double rho_c = 3.0*Mc / (4.0*pi*P3(R));
	
	double beta_j = beta(Gj); 

	double beta_c = beta(g*Gj);
	double beta_rel = (beta_j - beta_c) / (1.0 - beta_j*beta_c);
	double G_rel = 1.0 / sqrt(1.0 - P2(beta_rel));

	double rho_j = Lj / (P2(Gj)*pi*P2(jetRadius(z, theta))*beta_j*P3(cLight));

	double P_c = rho_j*P2(G_rel*beta_rel*cLight); //la presion de la nube es la jet ram pressure en el sistema del blob; hj=1

	double hc = 1.0 + adiaCoeff*P_c / ((adiaCoeff - 1.0)*rho_c*cLight2);
	
	double cs = sqrt(adiaCoeff*P_c / (hc*rho_c));
		
	double beta_cs = cs / cLight;
	if (beta_cs > 1.0) {
		std::cout << "error";
	}
	return cs;
}



void gammaC(State& st, Vector& Gc, Vector& Rc, Vector& tobs)
{
	//static const double Lj = GlobalConfig.get<double>("Lj");
	static const double theta = GlobalConfig.get<double>("openingAngle");
	static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double g = GlobalConfig.get<double>("g");
	static const double z_peak = GlobalConfig.get<double>("z_peak")*pc;
	static const double blobRadius = GlobalConfig.get<double>("blobRadius")
		                             *(z_peak*theta);

	//static const double Mc = GlobalConfig.get<double>("Mc")*solarMass;
	//static const double inc = GlobalConfig.get<double>("inc")*pi / 180;  //degree


	/*std::string archive = "L44_i60_short.txt";
	
	std::ifstream fileReaded;
	fileReaded.open(archive.c_str(), std::ios::in);
	double tobs_yr, t, z, g, Rc_Rj, rho, Q;
	*/

	
	/*double mu = cos(inc);

	const double z0 = st.electron.ps[DIM_R].first();
	const double zMax = st.electron.ps[DIM_R].last();
	const int N_R = st.electron.ps[DIM_R].size() - 1;

	double z_int = pow((zMax / z0), (1.0 / N_R));
		
	double R0 = stagnationPoint(z0);//aca dejo el R0, no el que se expande
	*/

	
	
	double beta_j = beta(Gj);

	st.electron.ps.iterate([&](const SpaceIterator& i) {

		//const double z = i.val(DIM_R);		
		int s = i.coord[DIM_R]; //posicion en la coordenada z

		//fileReaded >> tobs_yr >> t >> z >> g >> Rc_Rj >> rho >> Q;

		//z = z*pc;
		//tobs[s] = tobs_yr*yr;
		//Rc[s] = Rc_Rj*(z*theta);
		Gc[s] = g*Gj;
		Rc[s] = blobRadius;
		tobs[s] = 0.0;

	}, { 0, -1 }); //fijo cualquier energia
	//fileReaded.close();


	/*	double y = z / z0;

		
		double g;
		
		double beta_c;
		double beta_i = 1.0e-10;

		if (s == 0) {
			Rc[0] = R0;
			Gc[0] = 1.0;
			g = Gc[0] / Gj;
			tobs[0] = 0.0;
			beta_c = beta_i;
		}
		else {

			//primero calculo Gc[s]
			double rc = Rc[s - 1];
			Rc[s] = rc;
			
			double D = Lj*P2(rc) / (4.0*P2(theta)*P3(Gj*cLight)*z0*Mc); 

			g = Gc[s - 1] / Gj; // 1.0e-2; comienzo la iteracion con Gc el que entra


			double beta_j = sqrt(1.0 - 1.0 / P2(Gj));
			//double beta_c = 1.0 - 1.0 / (2.0*P2(Gc[s - 1]));

			double z_i = i.its[DIM_R].peek(-1);
			double dy = (z - z_i) / z0;
				
			
			if (s == 1) {
				beta_c = beta_i; // 1.0 - 1.0 / (2.0*(Gc(i_z - 1))**2);
			}
			else {
				beta_c = sqrt(1.0 - 1.0 / P2(Gc[s - 1]));
			}
			
			double dt = z0*dy / (beta_c*cLight);

			if (g < 1.0) {

				//double dgdy = P2(1.0 / g - g)* D / P2(y);  //*Gj estaba mal!
				//double dGdy = dgdy*Gj;						
				//Gc[s] = Gc[s-1]+dGdy*dy;

				double cte = Lj*P2(rc) / (Mc*cLight2 * P2(z*theta));
				double dGdt = cte*P2(1 - beta_c / beta_j) * beta_c*P2(Gc[s - 1]);
				Gc[s] = Gc[s - 1] + dGdt*dt;
			
			}
			else {
				Gc[s] = Gj;
			}

			//Calculo Rc[s]
			
			beta_c = beta(Gc[s]);
			dt = z0*dy / (beta_c*cLight);

			double beta_rel = (beta_j - beta_c) / (1.0 - beta_j*beta_c);
			double G_rel = 1.0 / sqrt(1.0 - P2(beta_rel));

			double Psn = P2(beta_rel*G_rel);
			double Plat = 1.0e-3*P2(beta_j*Gj);

			//double dt = z*(z_int - 1.0) / (beta_c*cLight);   // agregue el beta_c;
			   

			if (Plat < Psn) {
				//double cs = soundC(R0, z_int);
				Rc[s] = Rc[s-1] + cs*dt / Gc[s];
			}
			else { 
				double z_i = i.its[DIM_R].peek(-1);

				Rc[s] = Rc[s - 1]*z/z_i;

			}


			double cte = (1.0 - beta_c*mu);
			tobs[s] = tobs[s - 1] + dt*cte;
		}*/ 

	//}, { 0, -1 }); //fijo cualquier energia

	
}




double Fe(double g, double y)
{
	return pow(g, 4.0) * (1.0 / P2(g) - P2(g)) / P2(y);
}


double eEmax(double z, double Gc, double B, double Reff)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double accEfficiency = 1.0;// GlobalConfig.get<double>("accEfficiency"); //CAMBIAR A 0.1

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
	//static const double Gj = GlobalConfig.get<double>("Gamma");
	static const double z_peak = GlobalConfig.get<double>("z_peak")*pc;

	st.magf.fill([&](const SpaceIterator& i) {
		//double z = i.val(DIM_R);
		int z_ix = i.coord[DIM_R];

		//double beta_c = sqrt(1.0 - 1.0 / P2(Gc[z_ix]));
		//double beta_j = sqrt(1.0 - 1.0 / P2(Gj));
		//double beta_rel = (beta_j - beta_c) / (1.0 - beta_j*beta_c);
		//double G_rel = 1.0 / sqrt(1.0 - P2(beta_rel));

		return computeMagField(z_peak, Gc[z_ix]);  
	});
}


