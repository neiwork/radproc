#include "distribution.h"

#include "messages.h"

#include "injection.h"
#include "losses.h"
#include "timeDistribution.h"
//#include <fmath\interpolation.h>
//#include <fmath\elimiGaussiana.h>
//#include <fmath\matrixInit.h>
#include <fmath\physics.h>


void distribution(Particle& p, State& st)
{


	show_message(msgStart, Module_electronDistribution);


	ParamSpaceValues Naux(p.ps);


	for (size_t t_i = 0; t_i < p.ps[2].size(); t_i++) {

		int t_position = t_i; // posicion en la dimension z 

		for (size_t z_i = 0; z_i < p.ps[1].size(); z_i++) {

			int z_position = z_i; // posicion en la dimension z 

			//genero el psv auxiliar para el siguiente iterate
			Naux.fill([&p](const SpaceIterator& i){
				return p.distribution.get(i); //copia del N  
			}, { -1, z_position, t_position }); //copio el iterador que sigue asi los aux se completan bien


			p.distribution.fill([&p, &st, &Naux](const SpaceIterator& i){

				double Emax = pow(10.0, p.logEmax)*1.6e-12;

				double E = i.par.E;
				double t = i.par.E;
				double r = i.par.E;


				double Eeff = effectiveE(E, Emax, t, p, st);
				double dist1(0.0), dist2(0.0), dist3(0.0);

				if (i.its[1].canPeek(-1))
				{
					//estos son los puntos donde Q=0, y las particulas vienen de ti-1
					if (i.its[2].canPeek(-1))
					{
						double dist = Naux.interpolate({ Eeff, r, i.its[2].peek(-1) });
						double ratioLosses = losses(Eeff, p, st) / losses(E, p, st);
						dist2 = dist*ratioLosses;
					}
					else  //(t_position == 0)
					{
						dist2 = 0.0;
					}

				}
				else //(z_position == 0)
				{
					dist1 = timeDistribution(E, r, t, p, st, Eeff);
				}


				Naux.set(i, dist1 + dist2);


				if (i.its[1].canPeek(-1)) //if (z_position != 0)
				{
					double delta_t = 0.0; //VER
					double delta_xk = i.its[1].peek(0) - i.its[1].peek(-1); // z[k] - z[k - 1];

					SpaceCoord coord = i.moved({ 0, -1, 0 }); //N(ri-1)

					double ni_1 = Naux.get(coord) / delta_xk;
					double ni = 0.0;

					if (i.its[1].canPeek(1))
					{
						double delta_xk1 = i.its[1].peek(1) - i.its[1].peek(0); // z[k+1] - z[k];
						ni = Naux.get(i) / delta_xk1;
					}

					dist3 = (ni_1 - ni)*(delta_t*cLight);
				}
				else
				{
					dist3 = Naux.get(i); //dejo el que estaba de los pasos 1 y 2
				}

				return dist3;// p.distribution.fill(dist3);

			}, { -1, z_position, t_position });

		}//for sobre r

	}//for sobre t

	show_message(msgEnd, Module_electronDistribution);

}



