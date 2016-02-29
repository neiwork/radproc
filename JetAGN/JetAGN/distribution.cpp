#include "distribution.h"

#include "injection.h"
#include "losses.h"
#include "timeDistribution.h"
#include <fmath\interpolation.h>
#include <fmath\elimiGaussiana.h>
#include <fmath\matrixInit.h>
#include <fmath\physics.h>


void distribution(Particle& p, State& st)
{

	ParamSpaceValues Naux(p.ps);

	for (size_t z_i = 0; z_i < p.ps[1].size(); z_i++) {
		
	//	double z = p.ps[1][z_i]; //valor de la variable, lo comento porque con el z_position lo puedo calcular
		int z_position = z_i; // posicion en la dimension z 

	//p.ps.iterate([&p, &st](const SpaceIterator& j){  //reemplace este iterate con el for de arriba
		//double z = j.par.R; //valor de la variable que recorro en este iterate
		//int z_position = j.coord[1];

/*******************bloque 1*************************/

		double Eeff = 0.0;

		timeDistribution(p, st, z_position, Eeff);  

/*******************bloque 2*************************/

		//genero el psv auxiliar para el siguiente iterate
		Naux.fill([&p](const SpaceIterator& i){
			return p.distribution.get(i); //copia del N  //aca se rompe
		}, { -1, z_position, -1 }); //fijo z //copio el iterador que sigue asi los aux se completan bien

		p.ps.iterate([&p, &st, &z_position, &Eeff, &Naux](const SpaceIterator& i){

			double N2;
			if (i.its[2].canPeek(-1))
			{
				//dist = N1(Eeff, t[i - 1]);  
				double dist = p.distribution.interpolate({ Eeff, p.ps[1][z_position], i.its[2].peek(-1) });
				double E = i.par.E;
				double ratioLosses = losses(Eeff, p, st) / losses(E, p, st);
				
				N2 = dist*ratioLosses;
			}
			else{
				N2 = 0.0; //ver que pasa aca, sería en t[0], creo que esta bien que sea 0, no cambi ael N
			}

			double  Ntot = Naux.get(i) + N2;
			/*Naux.get(i)      = N1
			  dist*ratioLosses = N2 == N2(E[l], t[i]) = N1(Eeff, t[i - 1])*losses(Eeff, p, st) / losses(E[l], p, st);
			  distribution(E[l], z[k], t[i]) = N1(E[l], z[k], t[i]) + N2(E[l], z[k], t[i]);*/

			p.distribution.set(i, Ntot);

		}, { -1, z_position, -1 }); //fijo z

	} //cierro el for sobre el r(i)

//}, { 0, -1, 0 }); //afuera itero solo sobre z reemplace este iterate con un for sobre r(i)


/*******************bloque 3*************************/

	//reemplazo el siguiente iterate con dos for
	
	//p.ps.iterate([&p, &st](const SpaceIterator& j){  //en este fijo el r[0] //VER el doble iterate

		//int E_position = j.coord[0];
		//int t_position = j.coord[2];

		/*double delta_t;

		if (j.its[2].canPeek(-1)){

			delta_t = (j.par.T - j.its[2].peek(-1)) / Gamma;//    delta_t = (t[i] - t[i - 1]) / Gamma; 
		}	//std::cout << " (last index: " << i.coord.dims[0] - 1 << ", value: " << i.its[0].peek(-1) << ")";
		else{ //// este seria el t[0], en ese caso construyo el delta_t hacia adelante
			delta_t = (-j.par.T + j.its[2].peek(1)) / Gamma;
		}*/ 


	//genero el psv auxiliar para el siguiente iterate
	p.ps.iterate([&p, &Naux](const SpaceIterator& i){

		Naux.set(i, p.distribution.get(i)); //copia del N  tal cual esta

	}); 

	//las interpolaciones en el siguiente paso las hago todas sobre Naux, asi no interpolo con lo que se va cambiando

	for (int t_position = 0; t_position < (int)p.ps[2].size(); t_position++) {

		double delta_t;

		if (t_position != 0){
			delta_t = (p.ps[2][t_position] - p.ps[2][t_position - 1]) / Gamma;//    delta_t = (t[i] - t[i - 1]) / Gamma; 
		}
		else{ //// este seria el t[0], en ese caso construyo el delta_t hacia adelante
			delta_t = (p.ps[2][t_position + 1] - p.ps[2][t_position]) / Gamma;
		}

		for (int E_position = 0; E_position < (int)p.ps[0].size(); E_position++) {

			p.ps.iterate([&p, &st, &delta_t, &E_position, &t_position, &Naux](const SpaceIterator& i){

				double dist1, dist2;
				// example: gets value at the previous position of the iterator of dimension 0
				// (also checks if we can -- i.e. current index > 0)
				if (i.its[1].canPeek(-1))
				{

					double delta_xk = i.par.R - i.its[1].peek(-1); //delta_xk = z[k] - z[k - 1];

					if (i.its[1].canPeek(-2)){

						double delta_xk_1 = i.its[1].peek(-1) - i.its[1].peek(-2);// z[k - 1] - z[k - 2];

						dist1 = Naux.interpolate({ p.ps[0][E_position], i.its[1].peek(-1), p.ps[2][t_position] })
							*(delta_t*cLight / delta_xk_1);  //Ver

						dist2 = Naux.get(i)*(delta_t*cLight / delta_xk); //VER si me alcanza y no necesito aux

						p.distribution.set(i, dist1 - dist2);

						//distribution(E[l], z[k], t[i]) = distribution(E[l], z[k - 1], t[i])*(delta_t*cLight / delta_xk - 1)
						//							- distribution(E[l], z[k], t[i])*(delta_t*cLight / delta_xk);;
					}
					else{ //si no hay dos atras en r, el delta_xk_1 = delta_xk
						double delta_xk_1 = delta_xk; //VER

						dist1 = Naux.interpolate({ p.ps[0][E_position], i.its[1].peek(-1), p.ps[2][t_position] })
							*(delta_t*cLight / delta_xk_1);  //Ver

						dist2 = Naux.get(i)*(delta_t*cLight / delta_xk); //VER si me alcanza y no necesito aux

						p.distribution.set(i, dist1 - dist2);
					}
				}
				else {// este seria el z[0], no puedo mirar una para atras
					dist1 = 0.0; //VER que pasa	aca
					p.distribution.set(i, dist1);
				}

			}, { E_position, -1, t_position });

		} //cierro el for sobre E
	} //cierro el for sobre t
	
		
	//}, { -1, 0, -1 }); //VER necesito fijar alguna coordenada??? , { 0, 0, 0 }); //afuera itero solo sobre t
	
}


