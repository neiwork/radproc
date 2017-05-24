#include "timeDistribution.h"



#include "losses.h"
#include "injection.h"

#include <fmath\interpolation.h>
#include <fmath\elimiGaussiana.h>
#include <fmath\matrixInit.h>
#include <fmath\physics.h>

#include <fparameters\Dimension.h>

#include <iostream>
#include <fstream>


double effectiveE(double Ee, double Emax, double t, double r, Particle& p, State& st, const SpaceCoord& i)
{
	double Eeffmin, Eeffmax, Eeff_int, sum_Eeff, dEeff, dtau_eff;

	Eeffmin = Ee;
	Eeffmax = Emax;
	double Eeff = Eeffmin;

	double nEeff = 100;

	Eeff_int = pow((Eeffmax / Eeffmin), 1.0/nEeff);  //0.00001
	sum_Eeff = 0.0;

	while ((sum_Eeff <= t) && (Eeff <= Eeffmax))	{// for (int k=0; k<=100000; ++k)	{

		dEeff = Eeff*(Eeff_int - 1.0);
		dtau_eff = dEeff / losses(Eeff, r, p, st, i);
		sum_Eeff = sum_Eeff + dtau_eff;

		//	if (sum_Eeff > times[j]) goto

		Eeff = Eeff*Eeff_int;

	}

	return Eeff;

}

double timeDistribution(double Ee, double r, double t, Particle& p, State& st, double Eeff, const SpaceCoord& i)
{ 
	const double magneticField{ st.magf.get(i) };
	double Emax = p.emax();// eEmax(r, magneticField);
	double tmin = st.electron.ps[2][0]; //en teoria esto no es necesario, porque solo llamo a etsa funcion en el tmin
	//chequear que t == tmin

	double EpMin = Ee;   //Ep = Eprima
	double EpMax = Eeff;
	double Ep = EpMin;

	int nEp = 50;
	double Ep_int = pow((EpMax / EpMin), 1.0 / nEp);  //0.001

	double dEp(0.0), inj(0.0);

	double sum_Ep = 0.0;

	for (int l=0; l < nEp; ++l)	{

		dEp = Ep*(Ep_int-1.0);

		if (Ep > Emax){
			inj = 0.0;
		}
		else{
			inj = p.injection.interpolate({ { DIM_E, Ep }, { DIM_R, r }, { DIM_T, tmin} }); //paso los valores E,r,t en donde quiero evaluar Q
			//interpolDoble(Ep, Time[j]-tau, Ee, Time, injection);  //chequear que ande!
		}

		sum_Ep = sum_Ep + inj*dEp;

		Ep = Ep*Ep_int;

	}

	double perdidas = losses(Ee, r, p, st, i);

	double total = sum_Ep / perdidas;
				
	return total;

}
