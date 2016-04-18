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
	Naux.fill([&Naux](const SpaceIterator& i){
		return 0.0;
	});

	for (int t_ix = 0; t_ix < p.ps[2].size(); t_ix++) {

		//t_ix; posicion en la dimension z 
		//z_ix;  posicion en la dimension z 

		for (int z_ix = 0; z_ix <= t_ix; z_ix++) { //z_ix < p.ps[1].size(); z_ix++) {
			

		//for (int t_ix = 0; t_ix < p.ps[2].size(); t_ix++) {
			//genero el psv auxiliar para el siguiente iterate
			


			p.ps.iterate([&p, &st, &Naux, &z_ix, &t_ix](const SpaceIterator& i){

				double Emax = pow(10.0, p.logEmax)*1.6e-12;

				double E = i.par.E;
				double r = i.par.R;
				double t = i.par.T;


				double Eeff = effectiveE(E, Emax, t, r, p, st);
				double dist1(0.0), dist2(0.0), dist3(0.0);

				if (z_ix == 0)
				{
					dist1 = timeDistribution(E, r, t, p, st, Eeff);
				}
				else //if (i.its[1].canPeek(-1)) //z_ix != 0) //
				{
					//estos son los puntos donde Q=0, y las particulas vienen de ti-1
					//if (i.its[2].canPeek(-1)) //si r!=0 -> t!=0 asi que no necesito preguntar
					//{
					double dist = Naux.interpolate({ Eeff, r, i.its[2].peek(-1) });
					double ratioLosses = losses(Eeff, r, p, st) / losses(E, r, p, st);
					dist2 = dist*ratioLosses;
					//}
					//		else  //(t_position == 0)
					//		{
					//			dist2 = 0.0;
					//		}

				}

				Naux.set(i, dist1 + dist2);

				double ni = dist1 + dist2; //Naux.get(i)

				double ni_1;

				if (i.its[1].canPeek(-1)) //if (z_ix != 0)
				{
					SpaceCoord coord = i.moved({ 0, -1, 0 }); //N(ri-1)
					ni_1 = Naux.get(coord);// / delta_xk; //VER por que no dan lo mismo Naux y p.distribution

					dist3 = (ni_1 - ni);// las fracciones delta_t*cLight/delta_x == 1 siempre *(delta_t*cLight);
				}
				else
				{
					dist3 = ni;// p.injection.get(i)*i.its[2].peek(0) - ni;  //Q*t-ni, el Q*t seria el ni_1
				}

				
					//double timeStep = i.its[2].peek(0);

					//double delta_t = timeStep;// /Gamma;
					//double delta_xk = i.its[1].peek(0) - i.its[1].peek(-1); // z[k] - z[k - 1];

					//double ni_1 = p.distribution.get(coord) / delta_xk;
			
				p.distribution.set(i,dist3);// p.distribution.fill(dist3);
				
				Naux.ps.iterate([&p, &Naux](const SpaceIterator& j){
					Naux.set(j, p.distribution.get(j)); //copia del N  
				});// , { -1, z_ix, t_ix }); //copio el iterador que sigue asi los aux se completan bien

			} , { -1, z_ix, t_ix });

		}//for sobre r

	}//for sobre t

	show_message(msgEnd, Module_electronDistribution);

}



