#include "timeDistribution.h"



#include "losses.h"
#include "injection.h"

#include <fmath\interpolation.h>
#include <fmath\elimiGaussiana.h>
#include <fmath\matrixInit.h>
#include <fmath\physics.h>
#include <iostream>
#include <fstream>


double effectiveE(double Ee, double Emax, double t, double r, Particle& p, State& st)
{
	double Eeffmin, Eeffmax, Eeff_int, sum_Eeff, dEeff, dtau_eff;

	Eeffmin = Ee;
	Eeffmax = Emax;
	double Eeff = Eeffmin;

	Eeff_int = pow((Eeffmax / Eeffmin), 0.01);  //0.00001
	sum_Eeff = 0.0;

	while ((sum_Eeff <= t) && (Eeff <= Eeffmax))	{// for (int k=0; k<=100000; ++k)	{

		dEeff = Eeff*(Eeff_int - 1.0);
		dtau_eff = dEeff / losses(Eeff, r, p, st);
		sum_Eeff = sum_Eeff + dtau_eff;

		//	if (sum_Eeff > times[j]) goto

		Eeff = Eeff*Eeff_int;

	}

	return Eeff;

}

double timeDistribution(double Ee, double r, double t, Particle& p, State& st, double Eeff)
{ 

	//aca declaro las variables del bucle
	//double Eeff este no lo declaro porque lo paso como argumento
	//double Eeffmin, Eeffmax, Eeff_int, sum_Eeff, dEeff, dtau_eff;
	double dtau, tau;
	double EpMin, EpMax, Ep, Ep_int, sum_Ep, dEp;
	double EppMin, EppMax, Epp, Epp_int, sum_Epp, dEpp;
	double perdidas, inj, dist, t_cool, t_carac;

	
////////////////CON ESTE IF CONTROLO SI LAS PERDIDAS SON RELEVANTES FRENTE A Tadv+Tdec//////////////////////
	perdidas = losses(Ee, r, p, st);

	t_cool = Ee/perdidas;

	//t_esc = escapeTime(Ee, p);

	//comento el siguiente if porque Tesc-> inf 
	//por lo que las perdidas siempre dominan
	//if (t_cool <= t){//  (t_cool <= t_esc)){ //| (t_cool <= Time[j]) )	{   //perdidas importantes

///////////////BUSCO EL Eeff/////////////////////////////////////////////
	//double Eeff = effectiveE(Ee, Emax, t, p, st);
///////////////////YA ENCONTRE EL Eeff/////////////////////////////////////
//lo calculo para todo t porque lo necesito independientemente del siguiente if

	if (t_cool <= t){

///////////////////AHORA CALCULO LA INTEGRAL "DOBLE" QUE DA N(E,t)/////////////////////////////////////
		EpMin = Ee;   //Ep = Eprima
		EpMax = Eeff;
		Ep    = EpMin;

		Ep_int = pow((EpMax/EpMin),0.01);  //0.001

		sum_Ep = 0.0;

		for (size_t l=0; l <= 10; ++l)	{ //100

			dEp = Ep*(Ep_int-1.0);
					
////aca calculo la integral del exponente con variable Epp = E'' //////////////

			EppMin = Ee;   //Ep = Eprima
			EppMax = Ep;
			Epp    = EppMin; 

			Epp_int = pow((EppMax/EppMin),0.01);

			sum_Epp  = 0.0;
			//sum_dexp = 0.0;

			if (EppMax > EppMin){

				for (size_t n=0; n <= 10; ++n)	{ //100

					dEpp = Epp*(Epp_int-1.0);

					dtau = dEpp/losses(Epp, r, p, st);
					//dexp = dtau/escapeTime(Epp,p); 

					sum_Epp  = sum_Epp + dtau;
					//sum_dexp = sum_dexp + dexp;

					Epp = Epp*Epp_int;
	
				}
			}

			tau = sum_Epp;
			/////////////////// ya calcule el tau ////////////////////////////

			double teval = 0.0;
			double tmin = st.electron.ps[2][0]; //el siguiente if es para no interpolar fuera del rango

			if (t - tau < tmin){ teval = tmin;    }
			else{             	 teval = t - tau; }

			double Emax = eEmax(r, magneticField);

			if (Ep > Emax){
				inj = 0.0;
			}
			else{
				inj = p.injection.interpolate({ { DIM_E, Ep }, { DIM_R, r }, { DIM_T, t - tau } }); //paso los valores E,r,t en donde quiero evaluar Q
				//interpolDoble(Ep, Time[j]-tau, Ee, Time, injection);  //chequear que ande!
			}

			

			dist = inj;  //*exp(-sum_dexp);lo comento por lo del Tesc //tau/escapeTime(Ep,particle));

			sum_Ep = dist*dEp+sum_Ep;

			Ep = Ep*Ep_int;

		}
		//particle.distribution[j*(ne+1)+i] = sum_Ep/perdidas;
		double total = sum_Ep / perdidas;
				
		dist = total;

	} //aca termina la parte del if en donde tcool < T //las perdidas son importantes
			
	else {  //t_cool > t_esc or t_cool > t  

		t_carac = min(t_cool, t); // min(t_cool, t_esc); //min(t_esc, Time[j]));

		inj = p.injection.interpolate({ { DIM_E, Ee }, { DIM_R, r }, { DIM_T, t } });//  p.injection.get(i);// injection[j*(ne + 1) + i];  VER
	
		dist = inj*t_carac; // t_esc;  //t_carac;  //exp(-tau*escapeTime(Ep,particle));

		//particle.distribution[j*(ne+1)+i] = dist;
	//	return dist;

	}

	return dist;
}
