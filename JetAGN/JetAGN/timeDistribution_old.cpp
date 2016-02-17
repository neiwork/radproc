#include "timeDistribution.h"

#include <omp.h>

#include "losses.h"
#include "injection.h"

#include <fmath\interpolation.h>
#include <fmath\elimiGaussiana.h>
#include <fmath\matrixInit.h>
#include <fmath\physics.h>
#include <iostream>
#include <fstream>

//void timeDistribution(Particle& particle, State& state)
void timeDistribution(Vector& Ee, Vector& Time)
{

	Vector Ee = particle.energyPoints;
	Vector Time = particle.timePoints;

	double Emax = pow(10.0, particle.logEmax);

	Vector injection = particle.injection;

	int ne = Ee.size()-1;    //el -1 es para que sea el nEnergyPoints y 
	int nt = Time.size()-1;

//aca declaro las variables del bucle
	double Eeffmin, Eeffmax, Eeff, Eeff_int, sum_Eeff, dEeff, dtau_eff;
	double dtau,dexp,tau;
	double EpMin, EpMax, Ep, Ep_int, sum_Ep, dEp;
	double EppMin, EppMax, Epp, Epp_int, sum_Epp, dEpp, sum_dexp;
	double perdidas,inj,dist,t_cool,t_esc,t_carac;
	int l,p,i;

	int j;


	for (j=0; j < (int) Time.size() ; ++j)	{   //VER QUE PASA CON EL ULTIMO!!!!!!!!!

		state.setTime(j);  

		for (i=0; i < Ee.size() ; ++i)	{     //ver como acomodo los limites
//			std::cout << "T" << j << "/" << Time.size() << " E" << i << "/" << Ee.size() << std::endl;

////////////////CON ESTE IF CONTROLO SI LAS PERDIDAS SON RELEVANTES FRENTE A Tadv+Tdec//////////////////////
			perdidas = losses(Ee[i], particle, state);

			t_cool = Ee[i]/perdidas;

			t_esc = escapeTime(Ee[i], particle);

			if ( (t_cool <= t_esc) ){ //| (t_cool <= Time[j]) )	{   //perdidas importantes

///////////////BUSCO EL Eeff/////////////////////////////////////////////
				Eeffmin = Ee[i];
				Eeffmax = Emax;
				Eeff    = Eeffmin;

				Eeff_int = pow((Eeffmax/Eeffmin),0.01);  //0.00001
				sum_Eeff  = 0.0;

				while ( (sum_Eeff <= Time[j]) && (Eeff <= Eeffmax) )	{// for (int k=0; k<=100000; ++k)	{

					dEeff = Eeff*(Eeff_int-1.0);
					dtau_eff = dEeff/losses(Eeff, particle, state);
					sum_Eeff = sum_Eeff + dtau_eff;

				//	if (sum_Eeff > times[j]) goto

					Eeff = Eeff*Eeff_int;

				}	
	///////////////////YA ENCONTRE EL Eeff/////////////////////////////////////

	///////////////////AHORA CALCULO LA INTEGRAL "DOBLE" QUE DA N(E,t)/////////////////////////////////////
				EpMin = Ee[i];   //Ep = Eprima
				EpMax = Eeff;
				Ep    = EpMin;

				Ep_int = pow((EpMax/EpMin),0.01);  //0.001

				sum_Ep = 0.0;

				for (l=0; l <= 100; ++l)	{ //1000

					dEp = Ep*(Ep_int-1.0);
					
///////////////////	aca calculo la integral del exponente con variable Epp = E'' //////////////

					EppMin = Ee[i];   //Ep = Eprima
					EppMax = Ep;
					Epp    = EppMin; 

					Epp_int = pow((EppMax/EppMin),0.01);

					sum_Epp  = 0.0;
					sum_dexp = 0.0;

					if (EppMax > EppMin){

						for (p=0; p <= 100; ++p)	{

							dEpp = Epp*(Epp_int-1.0);

							dtau = dEpp/losses(Epp, particle, state);
							dexp = dtau/escapeTime(Epp,particle); 

							sum_Epp  = sum_Epp + dtau;
							sum_dexp = sum_dexp + dexp;

							Epp = Epp*Epp_int;
	
						}
					}

					tau = sum_Epp;
	/////////////////// ya calcule el tau ////////////////////////////
					inj = interpolDoble(Ep, Time[j]-tau, Ee, Time, injection);
					dist = inj*exp(-sum_dexp);  //tau/escapeTime(Ep,particle));

					sum_Ep = dist*dEp+sum_Ep;

					Ep = Ep*Ep_int;

				}
				particle.distribution[j*(ne+1)+i] = sum_Ep/perdidas;

			} //aca termina la parte del if en donde las perdidas son importantes

			else {  //t_cool > t_esc or t_cool > t

				t_carac = min(t_cool,t_esc); //min(t_esc, Time[j]));

				inj = injection[j*(ne+1)+i];
	
				dist = inj*t_esc;  //t_carac;  //exp(-tau*escapeTime(Ep,particle));

				particle.distribution[j*(ne+1)+i] = dist;

			}
		} //cierra el de E	
	} //cierra el de t
}
