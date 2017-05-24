#include "distribution.h"

#include "messages.h"

#include "injection.h"
#include "losses.h"
#include "timeDistribution.h"

#include <iostream>
#include <fmath\physics.h>
#include <fparameters\Dimension.h>
#include <fparameters\SpaceIterator.h>
#include <fparameters\parameters.h>
#include <boost/property_tree/ptree.hpp>

void distribution(Particle& p, State& st)
{
	static const double Gamma = GlobalConfig.get<double>("Gamma");

	show_message(msgStart, Module_electronDistribution);

	const ParamSpace& ps{ p.ps };
	
	const double rmin = ps[DIM_R].first();
	const double rmax = ps[DIM_R].last();
	const int nR = ps[DIM_R].size()-1;

	ParamSpaceValues N2(p.ps);
	N2.fill([&](const SpaceIterator& i){
		return 0.0;
	});

	ParamSpaceValues N12(p.ps);
	N12.fill([&](const SpaceIterator& i){
		return 0.0;
	});

	static const std::string injector = GlobalConfig.get<std::string>("injector");

	bool multiple = (injector == "multiple");
	bool single = (injector == "single");
	int superior;

	superior = p.ps[1].size(); //este es el caso de multiple


	for (int t_ix = 0; t_ix < p.ps[2].size(); t_ix++) {

		if (t_ix % 5 == 0) {
			std::cout << "time: " << t_ix << std::endl;
		}
		 
		//for (int z_ix = 0; z_ix < p.ps[1].size(); z_ix++) { //emisores para todo z
		//for (int z_ix = 0; z_ix <= t_ix; z_ix++) {	//unico emisor en z=0

		if (single) { superior = t_ix + 1; }
		//el +1 es para que el for incluya el t_ix

		for (int z_ix = 0; z_ix < superior; z_ix++) { 
		

			N2.ps.iterate([&](const SpaceIterator& i){
				N2.set(i, p.distribution.get(i));  //copia del N  
			});      //ver si quizas lo puedo copiar en t-1 que es donde lo necesito

			p.ps.iterate([&](const SpaceIterator& i){
				const double E = i.val(DIM_E);
				const double r = i.val(DIM_R);
				const double t = i.val(DIM_T);
				const double magf = st.magf.get(i);
				//
				// equivale a:
				//	i.its[0].val();
				//	i.its[0].peek(DIM_E);
				//
				// equivalia a:
				//	 ex i.par.T;


				double Emax = eEmax(r, magf);
				
				double tp = t / Gamma; //time in the FF

				if (E < 10.0*Emax){
					double Eeff = effectiveE(E, Emax, tp, r, p, st, i);
					double dist1(0.0), dist2(0.0);

					dist1 = timeDistribution(E, r, tp, p, st, Eeff, i);

					if (t_ix != 0)
					{	//estos son los puntos donde Q=0, y las particulas vienen de ti-1
						//if (i.its[2].canPeek(-1)) 

						double dist = N2.interpolate({ { DIM_E, Eeff }, { DIM_R, r }, { DIM_T, i.its[2].peek(-1) } });
						double ratioLosses = losses(Eeff, r, p, st, i) / losses(E, r, p, st, i);
						dist2 = dist*ratioLosses;
					}

					N12.set(i, dist1 + dist2); //lo cargo en N12 mientras interpolo de N2
					//N12.set(i, dist1); //lo cargo en N12 mientras interpolo de N2
				}
				else
				{
					N12.set(i, 0.0);  //si E>Emax entonces anulo la N
				}

			}, { -1, z_ix, t_ix });


			p.ps.iterate([&](const SpaceIterator& i){

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

				//p.distribution.set(i, ni); // prueba
				
			} , { -1, z_ix, t_ix });

		}//for sobre r

	}//for sobre t

	show_message(msgEnd, Module_electronDistribution);

}





