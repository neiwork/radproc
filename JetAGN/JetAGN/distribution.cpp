#include "distribution.h"

#include "messages.h"

#include "injection.h"
#include "losses.h"
#include "timeDistribution.h"
#include <fmath\physics.h>


void distribution(Particle& p, State& st)
{


	show_message(msgStart, Module_electronDistribution);


	ParamSpaceValues Naux(p.ps);


	for (int t_ix = 0; t_ix < p.ps[2].size(); t_ix++) {

		//t_ix; posicion en la dimension z 
		//z_ix;  posicion en la dimension z 

		for (int z_ix = 0; z_ix < p.ps[1].size(); z_ix++) {
			

			//genero el psv auxiliar para el siguiente iterate
			Naux.fill([&p](const SpaceIterator& i){
				return p.distribution.get(i); //copia del N  
			}, { -1, z_ix, t_ix }); //copio el iterador que sigue asi los aux se completan bien


			p.ps.iterate([&p, &st, &Naux](const SpaceIterator& i){  

				double Emax = pow(10.0, p.logEmax)*1.6e-12;

				double E = i.par.E;
				double r = i.par.R;
				double t = i.par.T;				


				double Eeff = effectiveE(E, Emax, t, p, st);
				double dist1(0.0), dist2(0.0), dist3(0.0);

				if (i.its[1].canPeek(-1)) //z_ix != 0) //
				{
					//estos son los puntos donde Q=0, y las particulas vienen de ti-1
					if (i.its[2].canPeek(-1)) //t_position != 0
					{
						double dist = Naux.interpolate({ Eeff, r, i.its[2].peek(-1) });
						double ratioLosses = losses(Eeff, p, st) / losses(E, p, st);
						dist2 = dist*ratioLosses;
					}
			//		else  //(t_position == 0)
			//		{
			//			dist2 = 0.0;
			//		}

				}
				else //(z_position == 0)
				{
					dist1 = timeDistribution(E, r, t, p, st, Eeff);
				}


				//Naux.set(i, dist1 + dist2);

				double Ni = dist1 + dist2;

				if (i.its[1].canPeek(-1)) //if (z_ix != 0)
				{
					double timeStep = (timeMax - timeMin) / nTimes; //VER porque ahora NO es el mismo para todos
					double delta_t = timeStep / Gamma;
					double delta_xk = i.its[1].peek(0) - i.its[1].peek(-1); // z[k] - z[k - 1];

					SpaceCoord coord = i.moved({ 0, -1, 0 }); //N(ri-1)

					double ni_1 = p.distribution.get(coord) / delta_xk;
				//	double ni_1_aux = Naux.get(coord) / delta_xk; VER por que no dan lo mismo Naux y p.distribution
					double ni = 0.0;

					if (i.its[1].canPeek(1))
					{
						double delta_xk1 = i.its[1].peek(1) - i.its[1].peek(0); // z[k+1] - z[k];
						double dum = Ni;// Naux.get(i);
						ni = dum / delta_xk1;
					}
					else
					{
						double r_int = pow((rmax / rmin), (1.0 / nR));
						double delta_xk1 = i.its[1].peek(0)*(r_int-1.0); // VER
						double dum = Ni;// Naux.get(i);
						ni =  dum/ delta_xk1;
						
					}
					if (ni != 0.0){
						double dum = 0.0;
					}
					dist3 = (ni_1 - ni)*(delta_t*cLight);
				}
				else
				{
					dist3 = Ni;// Naux.get(i); //dejo el que estaba de los pasos 1 y 2
				}

				p.distribution.set(i,dist3);// p.distribution.fill(dist3);

			}, { -1, z_ix, t_ix });

		}//for sobre r

	}//for sobre t

	show_message(msgEnd, Module_electronDistribution);

}



