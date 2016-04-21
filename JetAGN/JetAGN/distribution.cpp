#include "distribution.h"

#include "messages.h"

#include "injection.h"
#include "losses.h"
#include "timeDistribution.h"
#include <fmath\physics.h>


void distribution(Particle& p, State& st)
{
	show_message(msgStart, Module_electronDistribution);
	
	ParamSpaceValues N2(p.ps);
	N2.fill([&N2](const SpaceIterator& i){
		return 0.0;
	});

	ParamSpaceValues N12(p.ps);
	N12.fill([&N12](const SpaceIterator& i){
		return 0.0;
	});


	for (int t_ix = 0; t_ix < p.ps[2].size(); t_ix++) {
		 
		for (int z_ix = 0; z_ix <= t_ix; z_ix++) { //z_ix < p.ps[1].size(); z_ix++) {

			N2.ps.iterate([&p, &N2](const SpaceIterator& i){
				N2.set(i, p.distribution.get(i));  //copia del N  
			});      //ver si quizas lo puedo copiar en t-1 que es donde lo necesito

			p.ps.iterate([&p, &st, &N2, &N12, &z_ix, &t_ix](const SpaceIterator& i){

				double Emax = pow(10.0, p.logEmax)*1.6e-12;

				double E = i.par.E;
				double r = i.par.R;
				double t = i.par.T;
				
				double Eeff = effectiveE(E, Emax, t, r, p, st);
				double dist1(0.0), dist2(0.0);

				if (z_ix == 0)
				{
					dist1 = timeDistribution(E, r, t, p, st, Eeff);
				}

				if (t_ix != 0)
				{	//estos son los puntos donde Q=0, y las particulas vienen de ti-1
					//if (i.its[2].canPeek(-1)) 

					double dist = N2.interpolate({ Eeff, r, i.its[2].peek(-1) });
					double ratioLosses = losses(Eeff, r, p, st) / losses(E, r, p, st);
					dist2 = dist*ratioLosses;
				}

				N12.set(i, dist1 + dist2); //lo cargo en N12 mientras interpolo de N2

			}, { -1, z_ix, t_ix });

			//VER: aca creo qeu tengo que copiar algo

			p.ps.iterate([&p, &st, &N12, &z_ix, &t_ix](const SpaceIterator& i){

				double ni = N12.get(i); //el ni es que que obtengo con N12

				double dist3;

				double ri = i.its[1].peek(0);
				double rip1;

				if (i.its[1].canPeek(-1)) //if (z_ix != 0)
				{
					SpaceCoord coord = i.moved({ 0, -1, 0 }); //N(ri-1)
					double ni_1 = p.distribution.get(coord); //este lo calculo con el p.dist porque ya esta en r-1

					
					double ri_1 = i.its[1].peek(-1);
					
					if (i.its[1].canPeek(1)) { rip1 = i.its[1].peek(+1); }
					else{
						double r_int = pow((rmax / rmin), (1.0 / nR));
						rip1 = ri*r_int; 
					}

					dist3 = ni*(1.0 - ri / rip1) + ni_1*(ri_1 / ri);
				}
				else //z_ix=0
				{
					rip1 = i.its[1].peek(+1);
					dist3 = ni*(1.0 - ri / rip1);  
				}

				p.distribution.set(i,dist3); // lleno p.distribution e interpolo en N12
				
			} , { -1, z_ix, t_ix });

		}//for sobre r

	}//for sobre t

	show_message(msgEnd, Module_electronDistribution);

}



